"""
Filter management module for DustyPY.
Handles the loading, caching, and synthetic photometry calculations for instrument bandpasses
using the internal STScI Vega spectrum for Zero Point calibration.
"""

import os
import glob
from typing import Dict, Tuple, Optional, List
import warnings

import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.io import ascii
from synphot import SpectralElement

from ..utils.physics import simpson_integrate


class Filter:
    """
    Wrapper optimisé pour les transmissions d'instruments avec calcul dynamique 
    du Point Zéro basé sur le spectre de Vega interne.
    """

    _instance_cache: Dict[str, "Filter"] = {}

    def __new__(cls, bandpass_name: str, *args, **kwargs):
        if bandpass_name in cls._instance_cache:
            return cls._instance_cache[bandpass_name]
        return super().__new__(cls)

    def __init__(self, bandpass_name: str, n_int: int = 10000):
        """
        Initialise le filtre et calcule les constantes (ZP, Normalisations).
        """
        if not hasattr(self, '_initialized'):
            self.name = bandpass_name
            self.path = self._find_filter_path(bandpass_name)

            # 1. Chargement de la transmission (synphot)
            self.bandpass = SpectralElement.from_file(self.path)
            self._pivot_wavelength = self.bandpass.pivot().to(u.um).value
            
            # 2. Définition de la grille d'intégration haute résolution
            w_min = np.min(self.bandpass.waveset.value)
            w_max = np.max(self.bandpass.waveset.value)
            w_int_aa = np.linspace(w_min, w_max, n_int) * u.AA
            self.w_int_um = w_int_aa.to(u.um).value
            self.throughput = self.bandpass(w_int_aa).value

            # 3. Pré-calcul des dénominateurs de normalisation
            # Photon (lambda * T) et Energie (T)
            self.norm_photon, _ = simpson_integrate(self.w_int_um, self.throughput * self.w_int_um)
            self.norm_energy, _ = simpson_integrate(self.w_int_um, self.throughput)

            # 4. Calcul automatique du Point Zéro (Vega)
            self._setup_vega_zp()
            
            self._initialized = True

        Filter._instance_cache[bandpass_name] = self

    @classmethod
    def get(cls, name: str, n_int: int = 10000) -> "Filter":
        return cls(name, n_int=n_int)

    def _find_filter_path(self, name: str) -> str:
        """Localise le fichier FITS du filtre dans data/filter/"""
        base_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'filter')
        pattern = os.path.join(base_path, '**', f"{name}.fits")
        files = glob.glob(pattern, recursive=True)
        if not files:
            raise FileNotFoundError(f"Filtre '{name}.fits' non trouvé dans {base_path}")
        return files[0]

    def _setup_vega_zp(self):
        """
        Calcule le flux de Point Zéro (Fzp) en intégrant le spectre de Vega 
        fourni par le STScI situé dans data/filter/vega/vega.dat.
        """
        vega_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)), 
            'data', 'filter', 'vega', 'vega.dat'
        )
        
        if not os.path.exists(vega_path):
            warnings.warn(f"Spectre de Vega non trouvé à {vega_path}. Magnitude Vega impossible.")
            self.zp_flux = None
            return

        # Chargement Vega: col 1 = AA, col 2 = erg/s/cm2/A
        vega_data = np.loadtxt(vega_path)
        w_vega_um = vega_data[:, 0] * 1e-4
        f_vega_wm2um = vega_data[:, 1] * 10.0 # Conversion erg/s/cm2/A -> W/m2/um

        # Détermination du mode (SVO : Gaia=Photon, Reste=Energie)
        use_photon = any(tag in self.name.lower() for tag in ['gaia3', 'gaia'])
        
        # Interpolation de Vega sur la grille du filtre
        f_vega_interp = np.interp(self.w_int_um, w_vega_um, f_vega_wm2um)

        if use_photon:
            num, _ = simpson_integrate(self.w_int_um, f_vega_interp * self.throughput * self.w_int_um)
            self.zp_flux = num / self.norm_photon
        else:
            num, _ = simpson_integrate(self.w_int_um, f_vega_interp * self.throughput)
            self.zp_flux = num / self.norm_energy

    def calculate_synthetic_flux(
        self, 
        model_wavelength: np.ndarray, 
        model_flux: np.ndarray,
        use_photon_count: bool = True
    ) -> Tuple[float, float, float]:
        """
        Calcule le flux synthétique moyen à travers le filtre.
        """
        f_interp = np.interp(self.w_int_um, model_wavelength, model_flux)

        num_phot, _ = simpson_integrate(self.w_int_um, f_interp * self.throughput * self.w_int_um)
        num_ener, _ = simpson_integrate(self.w_int_um, f_interp * self.throughput)

        # Longueur d'onde effective (source-dependant)
        eff_wave = num_phot / num_ener if num_ener != 0 else self._pivot_wavelength

        # Flux moyen normalisé selon le mode
        if use_photon_count:
            f_eff = num_phot / self.norm_photon if self.norm_photon != 0 else 0.0
        else:
            f_eff = num_ener / self.norm_energy if self.norm_energy != 0 else 0.0

        return f_eff, 0.0, eff_wave

    def get_magnitude(self, flux_avg: float) -> float:
        """Calcule la magnitude Vega : -2.5 * log10(F_star / F_zp)"""
        if self.zp_flux is None:
            raise ValueError(f"Le Point Zéro de Vega n'est pas initialisé pour {self.name}.")
        return -2.5 * np.log10(flux_avg / self.zp_flux)

    @staticmethod
    def get_bandpass_data(filter_names: List[str]) -> Dict[str, Dict[str, np.ndarray]]:
        """
        Retrieves wavelength and throughput data for a list of filter names.

        Args:
        filter_names (list): List of filter names to retrieve.

        Returns:
        Dict[str, Dict[str, np.ndarray]]: Dictionary mapping filter names to
        dictionaries containing 'wavelength' (microns) and 'throughput' arrays.
        """
        results = {}

        for name in filter_names:
            # Récupère l'instance du filtre (utilise le cache interne)
            f = Filter.get(name)
            # On utilise les données pré-calculées lors de l'initialisation

            results[name] = {
            'wavelength': f.w_int_um,
            'throughput': f.throughput
            }
        return results

    @staticmethod
    def batch_compute(
        wavelength: np.ndarray, 
        flux: np.ndarray, 
        filter_names: List[str], 
        n_int: int = 10000, 
        wavelength_definition: str = "effective"
    ) -> Dict[str, Dict]:
        """
        Méthode de haut niveau pour calculer la photométrie sur plusieurs filtres.
        """
        results = {}
        for name in filter_names:
            f = Filter.get(name, n_int=n_int)
            
            # Logique hybride SVO
            use_photon = any(tag in name.lower() for tag in ['gaia3', 'gaia'])
            
            val, err, eff_wave = f.calculate_synthetic_flux(wavelength, flux, use_photon_count=use_photon)
            mag = f.get_magnitude(val)

            wave = eff_wave if wavelength_definition == "effective" else f._pivot_wavelength
            
            results[name] = {
                'wavelength': wave,
                'flux': val,
                'mag_vega': mag
            }
        return results