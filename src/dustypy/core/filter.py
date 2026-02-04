"""
Optimised filter class for DustyPY photometric calculations.
Manages loading, caching, metadata (detector type), and synthetic photometry
calculations with dynamic Vega calibration.
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
from astropy.units import UnitsWarning
from synphot import SpectralElement

# Importation des utilitaires DustyPY
from ..utils.physics import simpson_integrate


class Filter:
    """
    Optimised filter class for DustyPY photometric calculations.
    Manages loading, caching, metadata (detector type), and synthetic photometry
    calculations with dynamic Vega calibration.
    """

    _instance_cache: Dict[str, "Filter"] = {}
    
    # Path to the JSON metadata cache
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
        Initialize the Filter instance, loading metadata and preparing for photometric calculations.

        Args:
            bandpass_name (str): Name of the filter to load.
            n_int (int): Number of points for the high-resolution integration grid.
        """
        if not hasattr(self, '_initialized'):
            self.name = bandpass_name
            self.path = self._find_filter_path(bandpass_name)
            
            # 1. Load existing metadata from the JSON cache
            meta = self._load_metadata()
            self.detector_type = meta.get('detector_type', 'energy')
            self._pivot_wavelength = meta.get('pivot_um')
            self.zp_flux = meta.get('zp_wm2um')

            # 2. Load transmission using synphot (always needed for integration)
            self.bandpass = SpectralElement.from_file(self.path)
            
            # 3. Create high-resolution integration grid
            w_min = np.min(self.bandpass.waveset.value)
            w_max = np.max(self.bandpass.waveset.value)
            w_int_aa = np.linspace(w_min, w_max, n_int) * u.AA
            self.w_int_um = w_int_aa.to(u.um).value
            self.throughput = self.bandpass(w_int_aa).value

            # 4. Normalizations (Pre-calculated to speed up batch_compute)
            self.norm_photon, _ = simpson_integrate(self.w_int_um, self.throughput * self.w_int_um)
            self.norm_energy, _ = simpson_integrate(self.w_int_um, self.throughput)

            # 5. Calculate missing values and update cache if necessary
            needs_save = False
            
            if self._pivot_wavelength is None:
                self._pivot_wavelength = self.bandpass.pivot().to(u.um).value
                needs_save = True
            
            if self.zp_flux is None:
                self._setup_vega_zp()
                needs_save = True
            
            if needs_save:
                self._save_metadata()
            
            self._initialized = True

        Filter._instance_cache[bandpass_name] = self

    @classmethod
    def get(cls, name: str, n_int: int = 10000) -> "Filter":
        return cls(name, n_int=n_int)

    # --- METADATA MANAGEMENT (JSON) ---

    def _load_metadata(self) -> dict:
        """Retrieve the specific filter entry from the JSON file."""
        if os.path.exists(self.METADATA_PATH):
            with open(self.METADATA_PATH, 'r') as f:
                return json.load(f).get(self.name, {})
        return {}

    def _save_metadata(self):
        """Save the calculated filter characteristics to the JSON file."""
        all_meta = {}
        if os.path.exists(self.METADATA_PATH):
            with open(self.METADATA_PATH, 'r') as f:
                all_meta = json.load(f)
        
        all_meta[self.name] = {
            'detector_type': self.detector_type,
            'pivot_um': float(self._pivot_wavelength),
            'zp_wm2um': float(self.zp_flux) if self.zp_flux else None
        }

        with open(self.METADATA_PATH, 'w') as f:
            json.dump(all_meta, f, indent=4)

    @staticmethod
    def set_detector_type(filter_name: str, detector_type: str):
        """
        Sets the detector type, recalculates the ZP, and updates the JSON cache.

        Args:
            filter_name (str): Name of the filter to update.
            detector_type (str): 'photon' or 'energy'.
        """
        if detector_type not in ['photon', 'energy']:
            raise ValueError("detector_type must be 'photon' or 'energy'.")

        # Retrieve or create the instance to force recalculation
        f = Filter.get(filter_name)
        f.detector_type = detector_type
        
        # The ZP depends on the detector type, so we recalculate it necessarily
        f._setup_vega_zp()
        f._save_metadata()

        print(f"Detector type updated for {filter_name} -> {detector_type} (ZP recalculated)")

    # --- PHOTOMETRIC CALCULATIONS ---

    def _setup_vega_zp(self):
        """
        Calculate the Vega Zero Point flux (Fzp) by integrating the reference spectrum 
        located in data/filter/vega/vega.dat.
        """
        vega_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)), 
            'data', 'filter', 'vega', 'vega.dat'
        )
        if not os.path.exists(vega_path):
            self.zp_flux = None
            return

        # Load Vega: col 1 = AA, col 2 = erg/s/cm2/A
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
        """
        Calculate the synthetic flux for the given model spectrum.
        
        Args:
            model_wavelength (np.ndarray): Wavelength array of the model spectrum in microns.
            model_flux (np.ndarray): Flux array of the model spectrum in W/m²/μm.
            use_photon_count (bool): Whether to use photon-counting or energy integration.
        Returns:
            Tuple[float, float, float]: (flux_avg, flux_avg_error, effective_wavelength)
        """
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
        """
        Calculate the Vega magnitude: -2.5 * log10(F_star / F_zp).
        
        Args:
            flux_avg (float): Average flux of the object in W/m²/μm.
        Returns:
            float: Vega magnitude.
        """
        if self.zp_flux is None:
            raise ValueError(f"ZP not initialized for {self.name}. Check vega.dat.")
        return -2.5 * np.log10(flux_avg / self.zp_flux)

    # --- STATIC AND UTILITY METHODS ---

    @staticmethod
    def print_available_filter():
        """
        Print a formatted table of all available filters with their metadata.
        """
        if not os.path.exists(Filter.METADATA_PATH):
            print("No metadata available. Initialize at least one filter.")
            return

        with open(Filter.METADATA_PATH, 'r') as f:
            meta_all = json.load(f)
        
        filter_names = sorted(meta_all.keys())
        
        header = f"{'Filter':<25} | {'Type':<10} | {'Pivot [μm]':<12} | {'ZP [W/m²/μm]':<15}"
        print(f"\n{header}")
        print("-" * len(header))
        
        for name in filter_names:
            m = meta_all[name]
            type_det = m.get('detector_type', 'energy')
            pivot = f"{m.get('pivot_um'):.4f}" if m.get('pivot_um') else "N/A"
            zp = f"{m.get('zp_wm2um'):.4e}" if m.get('zp_wm2um') else "N/A"
            print(f"{name:<25} | {type_det:<10} | {pivot:<12} | {zp:<15}")
        print("-" * len(header) + "\n")

    @staticmethod
    def batch_compute(wavelength, flux, filter_names, n_int=10000, wavelength_definition="effective"):
        """
        Compute photometry for a list of filters respecting their detector_type.
        
        Args:
            wavelength (np.ndarray): Wavelength array of the model spectrum in microns.
            flux (np.ndarray): Flux array of the model spectrum in W/m²/μm.
            filter_names (List[str]): List of filter names to compute.
            n_int (int): Number of points for the high-resolution integration grid.
            wavelength_definition (str): "effective" or "pivot" for returned wavelength.
        Returns:
            Dict[str, Dict[str, float]]: Dictionary with filter names as keys and
            values containing 'wavelength', 'flux', 'error', and 'mag_vega'.
        """
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
        """
        Add a new filter and record its detector type in the manifest.
        
        Args:
            ascii_path (str): Path to the ASCII file containing wavelength and throughput.
            filter_name (str): Name to assign to the new filter.
            detector_type (str): 'photon' or 'energy'.
        Returns:
            Filter: The newly created Filter instance.
        """
        base_data_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'filter', 'custom')
        os.makedirs(base_data_path, exist_ok=True)
        fits_path = os.path.join(base_data_path, f"{filter_name}.fits")

        # Conversion ASCII -> FITS
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UnitsWarning)
            data = ascii.read(ascii_path)
            table = Table(data, names=('WAVELENGTH', 'THROUGHPUT'), units=('ANGSTROMS', 'TRANSMISSION'))
            table.write(fits_path, format='fits', overwrite=True)
            
            # Enregistrement initial (déclenchera le calcul du pivot/ZP à la première instanciation)
            cls.set_detector_type(filter_name, detector_type)

        print(f"Filter {filter_name} successfully added from {ascii_path}.")
        
        return cls(filter_name)

    def _find_filter_path(self, name: str) -> str:
        """
        Locate the physical FITS file.
        
        Args:
            name (str): Name of the filter.
        Returns:
            str: Path to the filter FITS file.
        Raises:
            FileNotFoundError: If the filter file is not found.
        """
        base_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'filter')
        pattern = os.path.join(base_path, '**', f"{name}.fits")
        files = glob.glob(pattern, recursive=True)
        if not files:
            raise FileNotFoundError(f"Filter '{name}' not found.")
        return files[0]
    
    @staticmethod
    def get_bandpass_data(filter_names: List[str]) -> Dict[str, Dict[str, np.ndarray]]:
        """
        Retrieve wavelength and throughput arrays for multiple filters.
        Args:
            filter_names (List[str]): List of filter names. 
        Returns:
            Dict[str, Dict[str, np.ndarray]]: Dictionary with filter names as keys and
            values containing 'wavelength' and 'throughput' arrays.
        """

        results = {}
        for name in filter_names:
            f = Filter.get(name)
            results[name] = {'wavelength': f.w_int_um, 'throughput': f.throughput}
        return results
    
    def get_pivot_wavelength(self) -> float:
        """
        Get the pivot wavelength of the filter in microns.

        Returns:
            float: Pivot wavelength in microns.
        """
        return self._pivot_wavelength
    
    def __str__(self):
        return f"Filter name = {self.name:<10} \ndetector_type = {self.detector_type:<10} \npivot_um = {self._pivot_wavelength:<10.4f}\nZP = {self.zp_flux:.4e} W/m²/μm {"":<10}"