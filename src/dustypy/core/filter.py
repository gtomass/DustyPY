"""
Module de gestion des filtres pour DustyPY.
Gère le chargement, le cache, les métadonnées (Type de détecteur) et les 
calculs de photométrie synthétique avec calibration Vega dynamique.
"""

import os
import json
import glob
from typing import Dict, Tuple, Optional, List
import warnings

import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.io import ascii
from synphot import SpectralElement

# Import de l'intégrateur local
from ..utils.physics import simpson_integrate


class Filter:
    """
    Wrapper optimisé pour les filtres instrumentaux.
    Gère les types de détecteurs (Photon/Energy) via un manifeste JSON.
    """

    _instance_cache: Dict[str, "Filter"] = {}
    
    # Chemin vers le fichier de métadonnées dans le package
    METADATA_PATH = os.path.join(
        os.path.dirname(os.path.dirname(__file__)), 
        'data', 'filter', 'metadata.json'
    )

    def __new__(cls, bandpass_name: str, *args, **kwargs):
        if bandpass_name in cls._instance_cache:
            return cls._instance_cache[bandpass_name]
        return super().__new__(cls)

    def __init__(self, bandpass_name: str, n_int: int = 10000):
        """
        Initialise le filtre, charge son type de détecteur et calcule le ZP Vega.
        """
        if not hasattr(self, '_initialized'):
            self.name = bandpass_name
            self.path = self._find_filter_path(bandpass_name)
            
            # 1. Charger le type de détecteur (Manifeste ou défaut)
            self.detector_type = self._load_detector_type()

            # 2. Charger la transmission via synphot
            self.bandpass = SpectralElement.from_file(self.path)
            self._pivot_wavelength = self.bandpass.pivot().to(u.um).value
            
            # 3. Créer la grille d'intégration
            w_min = np.min(self.bandpass.waveset.value)
            w_max = np.max(self.bandpass.waveset.value)
            w_int_aa = np.linspace(w_min, w_max, n_int) * u.AA
            self.w_int_um = w_int_aa.to(u.um).value
            self.throughput = self.bandpass(w_int_aa).value

            # 4. Normalisations (Pré-calculées pour accélérer batch_compute)
            self.norm_photon, _ = simpson_integrate(self.w_int_um, self.throughput * self.w_int_um)
            self.norm_energy, _ = simpson_integrate(self.w_int_um, self.throughput)

            # 5. Calculer le Point Zéro Vega dynamique
            self._setup_vega_zp()
            
            self._initialized = True

        Filter._instance_cache[bandpass_name] = self

    @classmethod
    def get(cls, name: str, n_int: int = 10000) -> "Filter":
        return cls(name, n_int=n_int)

    # --- GESTION DES MÉTADONNÉES (JSON) ---

    def _load_detector_type(self) -> str:
        """Récupère le type de détecteur dans le fichier JSON."""
        if os.path.exists(self.METADATA_PATH):
            with open(self.METADATA_PATH, 'r') as f:
                meta = json.load(f)
                return meta.get(self.name, {}).get('detector_type', 'energy')
        return 'energy'

    @staticmethod
    def set_detector_type(filter_name: str, detector_type: str):
        """
        Définit et sauvegarde le type de détecteur pour un filtre.
        Types valides : 'photon' ou 'energy'.
        """
        if detector_type not in ['photon', 'energy']:
            raise ValueError("detector_type doit être 'photon' ou 'energy'.")

        # Charger le manifeste existant
        meta = {}
        if os.path.exists(Filter.METADATA_PATH):
            with open(Filter.METADATA_PATH, 'r') as f:
                meta = json.load(f)

        # Mettre à jour
        if filter_name not in meta:
            meta[filter_name] = {}
        meta[filter_name]['detector_type'] = detector_type

        # Sauvegarder
        with open(Filter.METADATA_PATH, 'w') as f:
            json.dump(meta, f, indent=4)
        
        # Si l'instance est déjà en cache, on la met à jour
        if filter_name in Filter._instance_cache:
            Filter._instance_cache[filter_name].detector_type = detector_type
            Filter._instance_cache[filter_name]._setup_vega_zp() # Recalcul du ZP

        print(f"Type de détecteur mis à jour pour {filter_name} -> {detector_type}")

    # --- CALCULS PHOTOMÉTRIQUES ---

    def _setup_vega_zp(self):
        """Calcule le flux de Point Zéro Vega (Fzp) selon le mode du détecteur."""
        vega_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)), 
            'data', 'filter', 'vega', 'vega.dat'
        )
        if not os.path.exists(vega_path):
            self.zp_flux = None
            return

        vega_data = np.loadtxt(vega_path)
        w_vega_um = vega_data[:, 0] * 1e-4
        f_vega_wm2um = vega_data[:, 1] * 10.0 # erg/s/cm2/A -> W/m2/um
        f_vega_interp = np.interp(self.w_int_um, w_vega_um, f_vega_wm2um)

        if self.detector_type == 'photon':
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
        """Intégration du modèle sur le filtre."""
        f_interp = np.interp(self.w_int_um, model_wavelength, model_flux)

        num_phot, num_phot_err = simpson_integrate(self.w_int_um, f_interp * self.throughput * self.w_int_um)
        num_ener, num_ener_err = simpson_integrate(self.w_int_um, f_interp * self.throughput)

        eff_wave = num_phot / num_ener if num_ener != 0 else self._pivot_wavelength

        if use_photon_count:
            f_eff = num_phot / self.norm_photon if self.norm_photon != 0 else 0.0
            f_eff_err = num_phot_err / self.norm_photon if self.norm_photon != 0 else 0.0
        else:
            f_eff = num_ener / self.norm_energy if self.norm_energy != 0 else 0.0
            f_eff_err = num_ener_err / self.norm_energy if self.norm_energy != 0 else 0.0

        return f_eff, f_eff_err, eff_wave
    
    def get_magnitude(self, flux_avg: float) -> float:
        """Calcule la magnitude Vega : -2.5 * log10(F_star / F_zp)"""
        if self.zp_flux is None:
            raise ValueError(f"ZP non initialisé pour {self.name}. Vérifiez vega.dat.")
        return -2.5 * np.log10(flux_avg / self.zp_flux)

    # --- MÉTHODES STATIQUES ET UTILITAIRES ---

    @staticmethod
    def batch_compute(wavelength, flux, filter_names, n_int=10000, wavelength_definition="effective"):
        """Calcule la photométrie pour une liste de filtres en respectant leur detector_type."""
        results = {}
        for name in filter_names:
            f = Filter.get(name, n_int=n_int)
            use_photon = (f.detector_type == 'photon')
            
            val, err, eff_wave = f.calculate_synthetic_flux(wavelength, flux, use_photon_count=use_photon)
            mag = f.get_magnitude(val)
            
            wave = eff_wave if wavelength_definition == "effective" else f._pivot_wavelength
            results[name] = {
                'wavelength': wave, 
                'flux': val, 
                'error': err,        
                'mag_vega': mag
            }
        return results

    @classmethod
    def add_custom_filter(cls, ascii_path: str, filter_name: str, detector_type: str = 'energy'):
        """Ajoute un nouveau filtre et enregistre son type de détecteur."""
        base_data_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'filter', 'custom')
        os.makedirs(base_data_path, exist_ok=True)
        fits_path = os.path.join(base_data_path, f"{filter_name}.fits")

        # Conversion ASCII -> FITS
        data = ascii.read(ascii_path)
        table = Table(data, names=('WAVELENGTH', 'THROUGHPUT'), units=('ANGSTROMS', 'TRANSMISSION'))
        table.write(fits_path, format='fits', overwrite=True)
        
        # Enregistrement du detector_type
        cls.set_detector_type(filter_name, detector_type)
        
        return cls(filter_name)

    def _find_filter_path(self, name: str) -> str:
        base_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'filter')
        pattern = os.path.join(base_path, '**', f"{name}.fits")
        files = glob.glob(pattern, recursive=True)
        if not files:
            raise FileNotFoundError(f"Filtre '{name}' non trouvé.")
        return files[0]
    
    @staticmethod
    def get_bandpass_data(filter_names: List[str]) -> Dict[str, Dict[str, np.ndarray]]:
        results = {}
        for name in filter_names:
            f = Filter.get(name)
            results[name] = {'wavelength': f.w_int_um, 'throughput': f.throughput}
        return results
    
    def get_pivot_wavelength(self) -> float:
        """Retourne la longueur d'onde pivot du filtre en microns."""
        return self._pivot_wavelength