"""
DustyPY Result module.
Handles the storage, physical de-normalization, and analysis of Dusty 
simulation outputs, including synthetic photometry generation.
"""

import os
from typing import Dict, List, Optional, Tuple, Any
from importlib import resources

from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt

from ..core.model import Model
from ..utils.physics import get_common_filters
from .filter import Filter


class Result:
    """
    Stores and de-normalizes the results of a Dusty simulation.

    This class converts dimensionless Dusty outputs into physical units 
    (SI) based on the model's distance and luminosity. It also manages 
    optimized synthetic photometry calculations using a shared filter cache.

    Attributes:
        wavelength (np.ndarray): Wavelength grid in microns.
        model (Model): The source model used for the simulation.
        flux_norm (np.ndarray): Dimensionless flux from Dusty's .stb file.
        scalars (Dict[str, float]): Scalar values parsed from Dusty's .out file.
        flux (np.ndarray): Physical flux density in W/m^2/um.
        photometry (Dict[str, Dict[str, float]]): Calculated synthetic photometry.
    """

    # Class-level cache to avoid re-mapping filters for the same dataset object
    _dataset_filter_cache: Dict[int, Dict[str, Filter]] = {}

    def __init__(
        self,
        wavelength: np.ndarray,
        flux_norm: np.ndarray,
        scalars: Dict[str, float],
        model: Model
    ):
        """
        Initializes the Result object and performs physical de-normalization.

        Args:
            wavelength (np.ndarray): Wavelength array in microns.
            flux_norm (np.ndarray): Normalized flux array from Dusty.
            scalars (Dict[str, float]): Dictionary of output scalar parameters.
            model (Model): The Model instance associated with this result.
        """
        self.wavelength = wavelength
        self.model = model
        self.flux_norm = flux_norm
        self.scalars = scalars

        # Convert dimensionless Dusty output to physical units (W/m2/um)
        self.flux = self._de_normalize()

        # Dictionary for photometry: { 'Filter_Name': {'flux': val, 'error': val, 'wavelength': val} }
        self.photometry: Dict[str, Dict[str, float]] = {}

    @property
    def true_radius(self) -> float:
        """
        Calculates the physical inner radius (r1) in meters.

        Returns:
            float: Physical radius in meters.
        """
        r1_cm = self.scalars.get('r1(cm)', 0.0)
        # Scale by the square root of the luminosity ratio as per Dusty logic
        return r1_cm * np.sqrt(self.model.total_luminosity / 1e4) * 1e-2

    @property
    def true_mass_loss(self) -> float:
        """
        Calculates the physical mass-loss rate in solar masses per year (M_sun/yr).

        Returns:
            float: Physical mass-loss rate.
        """
        mdot_dusty = self.scalars.get('Mdot', 0.0)
        l_tot = self.model.total_luminosity
        rgd = self.model.gas_to_dust_ratio
        rho = self.model.dust_density

        # Scaling law: Mdot ~ L^0.75 * sqrt(RGD * rho)
        return mdot_dusty * (l_tot / 1e4)**0.75 * np.sqrt((rgd / 200.0) * (rho / 3.0))

    def _de_normalize(self) -> np.ndarray:
        """
        Converts the normalized Dusty flux to physical units (W/m^2/um).

        Returns:
            np.ndarray: Physical flux density array.
        """
        # Physical constants
        l_sun_watts = 3.828e26      # Solar luminosity in Watts
        pc_to_m = 3.0857e16         # Parsec to meters

        l_total = self.model.total_luminosity  # Total system L in L_sun
        d_m = self.model.distance * pc_to_m    # Distance in meters

        # 1. Calculate the physical bolometric flux at distance D (W/m^2)
        f_bol = (l_total * l_sun_watts) / (4 * np.pi * d_m**2)

        # 2. Convert: (dimensionless lambda*F_lambda) * f_bol / lambda
        # Dusty outputs lambda*F_lambda. We multiply by bolometric scale 
        # and divide by lambda to isolate F_lambda.
        return (self.flux_norm * f_bol) / self.wavelength

    def compute_photometry(self, filter_list: List[str]):
        """
        Calculates synthetic photometry for a specific list of filters.

        Results are stored in the self.photometry dictionary.

        Args:
            filter_list (List[str]): List of filter names to compute.
        """
        for fname in filter_list:
            if fname in self.photometry:
                continue

            try:
                # Fast retrieval via Filter singleton cache
                f_obj = Filter.get(fname)

                # High-speed calculation (approx 0.1ms per call)
                f_val, f_err = f_obj.calculate_synthetic_flux(
                    self.wavelength, 
                    self.flux
                )

                self.photometry[fname] = {
                    'wavelength': f_obj.get_pivot_wavelength,
                    'flux': f_val,
                    'error': f_err
                }
            except Exception as e:
                print(f"Error computing photometry for '{fname}': {e}")

    def compute_all_photometry(self, dataset: Any):
        """
        Automatically computes photometry for all filters present in a Dataset.

        Uses an internal cache to map VizieR filter names to local 
        bandpass files efficiently.

        Args:
            dataset (Dataset): The observational dataset containing filter info.
        """
        if dataset is None:
            return

        dataset_id = id(dataset)

        # 1. Initialize cache for this dataset if not already present
        if dataset_id not in Result._dataset_filter_cache:
            viz_filters = None
            if hasattr(dataset, 'table') and dataset.table is not None:
                if 'sed_filter' in dataset.table.colnames:
                    viz_filters = np.unique(dataset.table['sed_filter'].data)
            
            if viz_filters is None and hasattr(dataset, 'filters'):
                viz_filters = np.unique(dataset.filters)

            if viz_filters is None or len(viz_filters) == 0:
                print("Warning: Dataset does not contain filter information.")
                return

            # Scan local filter directory
            filter_dir = resources.files('dustypy').joinpath('data/filter/comp/nonhst')
            local_files = [f for f in os.listdir(str(filter_dir)) if f.endswith('.fits')]
            
            # Map observational names to package filenames
            mapping = get_common_filters(viz_filters, local_files)

            # Pre-instantiate Filter objects
            cached_filters = {
                viz_name: Filter(local_name) 
                for viz_name, local_name in mapping.items()
            }
            Result._dataset_filter_cache[dataset_id] = cached_filters

        # 2. Perform fast mathematical integration
        active_filters = Result._dataset_filter_cache[dataset_id]

        for viz_name, filt_obj in active_filters.items():
            try:
                syn_flux, syn_err = filt_obj.calculate_synthetic_flux(
                    self.wavelength, 
                    self.flux
                )

                self.photometry[viz_name] = {
                    'wavelength': filt_obj.get_pivot_wavelength,
                    'flux': syn_flux,
                    'error': syn_err
                }
            except Exception as e:
                print(f"Failed to compute photometry for {filt_obj.name}: {e}")

    def plot(self, ax: Optional[plt.Axes] = None, photometry: bool = True, filters: List[str] = None, filters_colors: List[str] = None, fill: bool = True, **kwargs) -> plt.Axes:
        """
        Plots the SED and synthetic photometry.

        Args:
            ax (Optional[plt.Axes]): Matplotlib axes. If None, creates a new figure.
            photometry (bool): Whether to plot synthetic photometry points.
            filters (List[str]): List of filter names to plot transmission curves for.
            filters_colors (List[str]): List of colors for the filter transmission curves.
            fill (bool): Whether to fill under the filter transmission curves.
            **kwargs: Arguments passed to ax.loglog.

        Returns:
            plt.Axes: The plot axes.
        """
        if ax is None:
            plt.figure(figsize=(10, 6))
            ax = plt.gca()

        ax.loglog(self.wavelength, self.flux, **kwargs)

        if photometry and self.photometry:
            # Plot synthetic points
            for _, data in self.photometry.items():
                ax.errorbar(
                    data['wavelength'], data['flux'], yerr=data['error'],
                    fmt='d', markersize=6, color='orange'
                )
            # Add proxy for legend
            ax.errorbar(
                [], [], fmt='d', markersize=6, 
                label='Synthetic Photometry', color='orange'
            )

        if filters is not None:
            # Plot filter transmission curves
            for fname in filters:
                try:
                    f_obj = Filter.get(fname)
                    color = filters_colors[filters.index(fname)] if filters_colors and len(filters_colors) > filters.index(fname) else 'gray'
                    ax.plot(
                        f_obj.w_int_um, f_obj.throughput * f_obj.w_int_um * np.max(self.flux) / np.max(f_obj.throughput),
                        label=f"Filter: {fname}", color=color
                    )
                    if fill:
                        ax.fill_between(
                            f_obj.w_int_um, 
                            0, f_obj.throughput * f_obj.w_int_um * np.max(self.flux) / np.max(f_obj.throughput),
                            alpha=0.2, color=color
                        )
                except Exception as e:
                    print(f"Failed to plot filter {fname}: {e}")

        ax.set_xlabel(r"Wavelength [$\mu$m]")
        ax.set_ylabel(r"Flux $F_\lambda$ [W/m²]")
        ax.legend()
        return ax

    def get_photometry(self) -> Dict[str, Dict[str, float]]:
        """Returns the calculated photometry dictionary."""
        return self.photometry

    def get_sed(self) -> Tuple[np.ndarray, np.ndarray]:
        """Returns the full SED (wavelength array, flux array)."""
        return self.wavelength, self.flux

    def get_outputs(self) -> Dict[str, float]:
        """Returns the scalar values extracted from Dusty's .out file."""
        return self.scalars
    
    def get_available_filters(self) -> Dict[str, str]:
        """Returns a dictionary of available synthetic photometry filters."""
        local_filenames = [
            f for f in os.listdir(
                str(resources.files('dustypy').joinpath('data/filter/comp/nonhst'))
            ) if f.endswith('.fits')
        ]

        available_bp = [
            os.path.splitext(f)[0].lower() for f in local_filenames
        ]

        return available_bp

    def __repr__(self) -> str:
        """Formatted summary of simulation results."""
        divider = "  " + "="*85 + "\n"
        label_w = 30
        
        info = "\n" + divider
        info += f"  {f'DustyPY Result: {self.model.name}':^85}\n"
        info += divider
        info += "  System Parameters:\n"
        info += (f"    {'Wavelength Range':<{label_w}}: "
                f"{self.wavelength.min():.3e} to {self.wavelength.max():.3e} μm\n")
        info += f"    {'Total Luminosity':<{label_w}}: {self.model.total_luminosity:.3e} L_sun\n"
        info += f"    {'Distance':<{label_w}}: {self.model.distance:.2f} pc\n"
        info += f"    {'True Radius (r1)':<{label_w}}: {self.true_radius:.3e} m\n"

        if self.scalars:
            keys_to_show = [
                ('tau0', 'Optical Depth (tau0)'), 
                ('T1(K)', 'Inner Dust Temperature'),
                ('Td(K)', 'Outer Dust Temperature'), 
                ('Mdot', 'Mass-Loss Rate')
            ]
            for key, label in keys_to_show:
                if key in self.scalars:
                    val = self.true_mass_loss if key == 'Mdot' else self.scalars[key]
                    info += f"    {label:<{label_w}}: {val:.3e}\n"

        if self.photometry:
            info += divider
            # Calculate column size based on longest filter name
            col_size = max(len(name) for name in self.photometry.keys())
            col_size = max(col_size, 15)
            
            info += (f"  {'Filter':<{col_size}} | {'Pivot [μm]':<{col_size}} | "
                    f"{'Flux [W/m²]':<{col_size}} | {'Error [W/m²]':<{col_size}}\n")
            info += "  " + "-"*85 + "\n"
            
            # Sort by wavelength
            sorted_filters = sorted(
                self.photometry.items(), 
                key=lambda x: x[1]['wavelength']
            )
            for fname, data in sorted_filters:
                info += (f"  {fname:<{col_size}} | {data['wavelength']:<{col_size}.3e} | "
                        f"{data['flux']:<{col_size}.3e} | {data['error']:<{col_size}.3e}\n")

        info += divider
        return info