"""
Physical model module for DustyPY.
Defines the Model class which represents the physical state of the system, 
including stellar parameters, dust properties, and geometry.
"""

from typing import List, Dict, Optional, Any, Tuple
import numpy as np

from ..utils.physics import log_space, validate_wavelength_spacing


class Model:
    """
    Physical representation of the astronomical system.

    This class stores the configuration for a Dusty simulation, including 
    central stars, dust envelope properties, and spatial geometry. It also 
    tracks system-wide parameters like distance and total luminosity for 
    SED de-normalization.

    Attributes:
        name (str): Unique name for the model (used for file generation).
        distance (float): Distance to the system in parsecs.
        spectral_shape (str): Dusty spectral shape identifier (e.g., 'black_body').
        spectral_file (Optional[str]): Path to a custom spectral file for Dusty.
        grid_path (Optional[str]): Path to an atmosphere grid FITS file.
        stars (List[Dict[str, float]]): List of stellar components.
            Each dict contains 'temperature', 'luminosity', and 'logg'.
        dust (Dict[str, Any]): Dictionary containing dust chemical and 
            physical properties.
        gas_to_dust_ratio (float): The Gas-to-Dust mass ratio (RGD).
        dust_density (float): Bulk density of the dust grains in g/cm3.
        geometry (Dict[str, Any]): Dictionary defining the envelope geometry.
        wavelengths (Optional[np.ndarray]): Custom wavelength grid for 
            the simulation.
        dataset (Optional[Any]): Associated observational Dataset object.
    """

    def __init__(
        self,
        name: str = "default_model",
        distance: float = 100.0,
        gas_to_dust_ratio: float = 200.0,
        dust_density: float = 3.0
    ):
        """
        Initializes a new Model instance.

        Args:
            name (str): Name of the model.
            distance (float): Distance in parsecs. Defaults to 100.0.
            gas_to_dust_ratio (float): RGD ratio. Defaults to 200.0.
            dust_density (float): Grain density in g/cm3. Defaults to 3.0.
        """
        self.name = name
        self.distance = distance
        self.spectral_shape = "black_body"
        self.spectral_file: Optional[str] = None
        self.grid_path: Optional[str] = None
        
        # List of stars with absolute luminosity (L_sun) and logg
        self.stars: List[Dict[str, float]] = []

        self.dust: Dict[str, Any] = {
            'opacity': 0.1,
            'temperature': 300.0,
            'sio_absorption': 10.0,
            'sublimation_temp': 1500.0,
            'properties_index': 'common_and_addl_grain_composite',
            'composition': {'sil-dl': 1.0},
            'size_distribution': {
                'type': 'MODIFIED_MRN', 
                'amin': 0.005, 
                'amax': 0.25, 
                'q': 3.5
            }
        }

        self.gas_to_dust_ratio = gas_to_dust_ratio
        self.dust_density = dust_density

        self.geometry: Dict[str, Any] = {
            'type': 'POWD',
            'thickness': 1000.0,
            'n_powers': 1,
            'power': 2.0
        }

        self.wavelengths: Optional[np.ndarray] = None
        self.dataset: Optional[Any] = None

    def add_star(
        self, 
        temperature: float, 
        luminosity: float = 1.0, 
        logg: float = 4.0
    ):
        """
        Adds a stellar component to the system.

        Args:
            temperature (float): Effective temperature in Kelvin.
            luminosity (float): Absolute luminosity in L_sun. Defaults to 1.0.
            logg (float): Surface gravity (log10). Defaults to 4.0.
        """
        self.stars.append({
            'temperature': temperature,
            'luminosity': luminosity,
            'logg': logg
        })

    def set_dust(
        self,
        opacity: Optional[float] = None,
        temperature: Optional[float] = None,
        composition: Optional[Dict[str, float]] = None,
        rgd: Optional[float] = None,
        density: Optional[float] = None,
        size_distribution: Optional[Dict[str, Any]] = None,
        sio_absorption: Optional[float] = None,
        sublimation_temp: Optional[float] = None,
        properties_index: Optional[str] = None,
        envelope_type: Optional[str] = None,
        thickness: Optional[float] = None,
        n_powers: Optional[int] = None,
        power: Optional[float] = None
    ):
        """
        Configures dust and geometry properties of the envelope.

        Args:
            opacity (Optional[float]): Optical depth at 0.55 microns (tau_V).
            temperature (Optional[float]): Dust temperature at the inner 
                boundary in Kelvin.
            composition (Optional[Dict[str, float]]): Grain species and 
                their relative abundances.
            rgd (Optional[float]): Gas-to-Dust mass ratio.
            density (Optional[float]): Bulk density of grains (g/cm3).
            size_distribution (Optional[Dict[str, Any]]): Grain size dist. parameters.
            sio_absorption (Optional[float]): Depth of SiO absorption (%).
            sublimation_temp (Optional[float]): Sublimation temperature (K).
            properties_index (Optional[str]): Dusty index for optical properties.
            envelope_type (Optional[str]): Geometry type (e.g., 'POWD', 'SHELL').
            thickness (Optional[float]): Relative thickness of the envelope (Y).
            n_powers (Optional[int]): Number of power-law components.
            power (Optional[float]): Power-law index for the density profile.
        """
        # Dust dictionary updates
        if opacity is not None:
            self.dust['opacity'] = opacity
        if temperature is not None:
            self.dust['temperature'] = temperature
        if composition is not None:
            self.dust['composition'] = composition
        if size_distribution is not None:
            self.dust['size_distribution'] = size_distribution
        if sio_absorption is not None:
            self.dust['sio_absorption'] = sio_absorption
        if sublimation_temp is not None:
            self.dust['sublimation_temp'] = sublimation_temp
        if properties_index is not None:
            self.dust['properties_index'] = properties_index
        
        # Geometry dictionary updates
        if envelope_type is not None:
            self.geometry['type'] = envelope_type
        if thickness is not None:
            self.geometry['thickness'] = thickness
        if n_powers is not None:
            self.geometry['n_powers'] = n_powers
        if power is not None:
            self.geometry['power'] = power
            
        # Global attributes
        if rgd is not None:
            self.gas_to_dust_ratio = rgd
        if density is not None:
            self.dust_density = density

    def set_distance(self, distance_pc: float):
        """
        Sets the system distance.

        Args:
            distance_pc (float): Distance in parsecs.
        """
        self.distance = distance_pc

    def set_dataset(self, dataset: Any):
        """
        Associates an observational Dataset with this model.

        Args:
            dataset (Dataset): The observational dataset object.
        """
        self.dataset = dataset

    @property
    def total_luminosity(self) -> float:
        """
        Calculates the total integrated absolute luminosity.

        Returns:
            float: Sum of stellar luminosities in solar units (L_sun).
        """
        if not self.stars:
            return 1.0
        return sum(star['luminosity'] for star in self.stars)

    def set_wavelength_grid(self, intervals: List[Tuple[float, float, int]]):
        """
        Defines a custom wavelength grid based on logarithmic intervals.

        Args:
            intervals (List[Tuple[float, float, int]]): A list of tuples, 
                each containing (start_wavelength, end_wavelength, num_points).
        """
        wl_points = []
        for start, end, num in intervals:
            wl_points.extend(log_space(start, end, num))
        
        # Sort and remove duplicates
        unique_wl = np.unique(np.sort(wl_points))
        
        # Validate grid spacing for Dusty compatibility
        if validate_wavelength_spacing(unique_wl):
            self.wavelengths = unique_wl

    def __repr__(self) -> str:
        """
        Returns a string representation of the Model.
        """
        return (
            f"<Model(name='{self.name}', "
            f"distance={self.distance}pc, "
            f"stars={len(self.stars)}, "
            f"L_tot={self.total_luminosity:.2e}L_sun)>"
        )