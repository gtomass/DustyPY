"""
Filter management module for DustyPY.
Handles the loading, caching, and synthetic photometry calculations for instrument bandpasses.
"""

import os
import glob
from typing import Dict, Tuple, Optional
import warnings

import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.io import ascii
from synphot import SpectralElement

from ..utils.physics import simpson_integrate


class Filter:
    """
    Optimized wrapper for instrument bandpasses with caching and pre-computation.

    This class provides high-performance synthetic flux calculations by 
    pre-calculating the filter's integration grid and normalization factor 
    upon initialization. It uses a singleton-like pattern to ensure each 
    unique filter is only loaded once into memory.

    Attributes:
        name (str): The name of the bandpass/filter.
        path (str): Absolute path to the filter's FITS file.
        bandpass (SpectralElement): The underlying synphot spectral element.
        w_int_um (np.ndarray): Pre-calculated integration wavelength grid in microns.
        throughput (np.ndarray): Filter transmission values on the integration grid.
        norm_factor (float): Pre-calculated normalization factor (denominator) 
            for photon-counting integration.
    """

    _instance_cache: Dict[str, "Filter"] = {}

    def __new__(cls, bandpass_name: str, *args, **kwargs):
        """
        Ensures only one instance of a specific filter exists in memory.

        Args:
            bandpass_name (str): The unique identifier/name of the filter.

        Returns:
            Filter: The cached or newly created Filter instance.
        """
        if bandpass_name in cls._instance_cache:
            return cls._instance_cache[bandpass_name]
        return super().__new__(cls)

    def __init__(self, bandpass_name: str, n_int: int = 5000):
        """
        Initializes the Filter instance and performs pre-calculations.

        Args:
            bandpass_name (str): Name of the bandpass (must match a .fits filename).
            n_int (int): Number of points for the integration grid. Defaults to 5000.
        """
        # Prevent re-initialization if retrieved from cache
        if hasattr(self, '_initialized'):
            return

        self.name = bandpass_name
        self.path = self._find_filter_path(bandpass_name)

        # 1. Load the filter file using synphot
        self.bandpass = SpectralElement.from_file(self.path)
        self._pivot_wavelength = self.bandpass.pivot().to(u.um).value

        # 2. Pre-calculate the fixed integration grid in microns
        w_min = np.min(self.bandpass.waveset.value)
        w_max = np.max(self.bandpass.waveset.value)
        w_int_aa = np.linspace(w_min, w_max, n_int) * u.AA

        self.w_int_um = w_int_aa.to(u.um).value
        self.throughput = self.bandpass(w_int_aa).value

        # 3. Pre-calculate the normalization denominator
        # Uses photon-counter logic: integral(Transmission * lambda * dlambda)
        self.norm_factor, _ = simpson_integrate(
            self.w_int_um, 
            self.throughput * self.w_int_um
        )

        self._initialized = True
        Filter._instance_cache[bandpass_name] = self

    @classmethod
    def get(cls, name: str) -> "Filter":
        """
        Utility method to access the filter cache.

        Args:
            name (str): Name of the filter to retrieve.

        Returns:
            Filter: The requested Filter instance.
        """
        return cls(name)

    def _find_filter_path(self, name: str) -> str:
        """
        Recursively searches for the filter's FITS file in the package data directory.

        Args:
            name (str): Name of the filter file (without extension).

        Returns:
            str: Absolute path to the found FITS file.

        Raises:
            FileNotFoundError: If no FITS file matches the provided name.
        """
        base_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)), 
            'data', 'filter'
        )
        pattern = os.path.join(base_path, '**', f"{name}.fits")
        files = glob.glob(pattern, recursive=True)
        
        if not files:
            raise FileNotFoundError(
                f"Filter file '{name}.fits' not found in {base_path}"
            )
        return files[0]

    @property
    def central_wavelength(self) -> float:
        """
        Returns the pivot wavelength of the filter.

        Returns:
            float: Pivot wavelength in microns.
        """
        return self._pivot_wavelength

    def calculate_synthetic_flux(
        self, 
        model_wavelength: np.ndarray, 
        model_flux: np.ndarray
    ) -> Tuple[float, float]:
        """
        Calculates the synthetic flux by integrating the model SED over the filter.

        This method interpolates the provided model onto the filter's pre-calculated 
        internal grid and performs a Simpson integration.

        Args:
            model_wavelength (np.ndarray): Wavelength array of the model (microns).
            model_flux (np.ndarray): Flux density array of the model (W/m2/um).

        Returns:
            Tuple[float, float]: A tuple containing (synthetic_flux, integration_error).
                Flux is returned in the same physical units as the input model_flux.
        """
        # Fast interpolation of the model onto the filter's high-res internal grid
        f_interp = np.interp(self.w_int_um, model_wavelength, model_flux)

        # Calculate numerator: integral(Flux * Transmission * lambda * dlambda)
        num, num_err = simpson_integrate(
            self.w_int_um, 
            f_interp * self.throughput * self.w_int_um
        )

        # The effective flux is the ratio of the integrals
        if self.norm_factor == 0:
            return 0.0, 0.0

        f_eff = num / self.norm_factor
        f_eff_err = num_err / self.norm_factor

        return f_eff, f_eff_err
    
    @staticmethod
    def batch_compute(wavelength: np.ndarray, flux: np.ndarray, filter_names: list):
        """
        Helper to compute photometry for multiple filters on arbitrary arrays.
        
        Args:
            wavelength (np.ndarray): Wavelength array of the model (microns).
            flux (np.ndarray): Flux density array of the model (W/m2/um).
            filter_names (list): List of filter names to compute photometry for.
        Returns:
            Dict[str, Dict]: Dictionary with filter names as keys and dictionaries
                             containing 'wavelength', 'flux', and 'error' as values.
        """
        results = {}
        for name in filter_names:
            f = Filter.get(name)
            val, err = f.calculate_synthetic_flux(wavelength, flux)
            results[name] = {'wavelength': f.central_wavelength, 'flux': val, 'error': err}
        return results
    
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
    
    @classmethod
    def add_custom_filter(cls, ascii_path: str, filter_name: str) -> "Filter":
        """
        Converts an ASCII filter file (e.g., from SVO) to FITS and adds it to the library.

        Args:
            ascii_path (str): Path to the source ASCII file.
            filter_name (str): Name to give the filter (e.g., 'my_custom_g').

        Returns:
            Filter: The newly created and loaded Filter instance.
        """
        # 1. Determine the internal destination path
        base_data_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)), 
            'data', 'filter', 'comp', 'nonhst'
        )
        os.makedirs(base_data_path, exist_ok=True)
        fits_path = os.path.join(base_data_path, f"{filter_name}.fits")

        # Check if the file already exists
        if os.path.exists(fits_path):
            print(f"Filter '{filter_name}' already exists in the library.")
            return cls(filter_name)

        # 2. Read ASCII and convert to FITS Table
        # Standard SVO files are usually space/tab separated with Wavelength/Throughput
        data = ascii.read(ascii_path)
        
        # We ensure names and units match synphot expectations
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=u.UnitsWarning)

            table = Table(
                data, 
                names=('WAVELENGTH', 'THROUGHPUT'), 
                units=('ANGSTROMS', 'TRANSMISSION')
            )
            
            # 3. Save to the library folder
            table.write(fits_path, format='fits', overwrite=True)
        print(f"Successfully added '{filter_name}' to DustyPY filter library.")

        # 4. Return the initialized instance
        return cls(filter_name)