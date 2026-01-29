"""
Dataset management module for DustyPY.
Handles the retrieval, aggregation, and physical conversion of observational SED data.
"""

import os
from typing import Optional, List, Dict, Union

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.units as u

from ..utils.physics import apply_extinction


class Dataset:
    """
    Handles observational data (wavelengths, fluxes, and errors).

    This class provides tools to store, manipulate, and convert photometric
    data points, including fetching them directly from the VizieR SED service.

    Attributes:
        wavelength (np.ndarray): Wavelength array in microns.
        flux (np.ndarray): Flux density array in physical units (default W/m2/um).
        flux_err (np.ndarray): Flux error array.
        filters (np.ndarray): Array of filter names/strings associated with data points.
        units (str): Current unit string for the flux (e.g., 'W/m2' or 'Jy').
        table (Optional[Table]): The underlying astropy Table if loaded from a catalog.
        filters_mapped (Dict[str, str]): Mapping of VizieR filter names to local bandpasses.
    """

    def __init__(
        self,
        wavelength: Optional[np.ndarray] = None,
        flux: Optional[np.ndarray] = None,
        flux_err: Optional[np.ndarray] = None,
        filters: Optional[np.ndarray] = None,
        units: str = "W/m2",
        table: Optional[Table] = None
    ):
        """
        Initializes the Dataset object.

        Args:
            wavelength (Optional[np.ndarray]): Initial wavelengths in microns.
            flux (Optional[np.ndarray]): Initial fluxes.
            flux_err (Optional[np.ndarray]): Initial flux errors.
            filters (Optional[np.ndarray]): Initial filter names.
            units (str): Units of the input flux. Defaults to "W/m2".
            table (Optional[Table]): Optional Astropy table for metadata.
        """
        self.wavelength = np.array(wavelength) if wavelength is not None else np.array([])
        self.flux = np.array(flux) if flux is not None else np.array([])
        self.flux_err = (np.array(flux_err) if flux_err is not None 
                         else np.zeros_like(self.flux))
        self.filters = (np.array(filters, dtype=str) if filters is not None 
                        else np.zeros_like(self.flux, dtype=str))
        self.units = units
        self.table = table
        self.filters_mapped: Dict[str, str] = {}

    def add_data(
        self,
        wavelength: Union[float, List[float], np.ndarray],
        flux: Union[float, List[float], np.ndarray],
        flux_err: Optional[Union[float, List[float], np.ndarray]] = None,
        filters: Optional[Union[str, List[str], np.ndarray]] = None
    ):
        """
        Appends new data points to the dataset.

        Args:
            wavelength: Wavelength(s) to add.
            flux: Flux value(s) to add.
            flux_err: Optional error value(s).
            filters: Optional filter name(s).
        """
        self.wavelength = np.concatenate([self.wavelength, np.atleast_1d(wavelength)])
        self.flux = np.concatenate([self.flux, np.atleast_1d(flux)])

        err = (np.atleast_1d(flux_err) if flux_err is not None 
               else np.zeros_like(np.atleast_1d(flux)))
        self.flux_err = np.concatenate([self.flux_err, err])

        filt = (np.atleast_1d(filters) if filters is not None 
                else np.zeros_like(np.atleast_1d(flux), dtype=str))
        self.filters = np.concatenate([self.filters, filt])

    def convert_to_jansky(self):
        """
        Converts flux and error from W/m2/um to Jansky (Jy).

        Conversion formula: F_nu [Jy] = 1e26 * F_lambda * (lambda^2 / c)
        where lambda is in meters and F_lambda is in W/m2/m.
        """
        if self.units == "W/m2":
            speed_of_light = 2.99792458e8  # m/s
            # Factor = 1e26 * (wl_um^2 * 1e-12) / (c * 1e-6)
            factor = 1e26 * (self.wavelength**2 * 1e-6) / speed_of_light
            
            self.flux *= factor
            self.flux_err *= factor
            self.units = "Jy"

    def convert_to_watt(self):
        """
        Converts flux and error from Jansky (Jy) to W/m2/um.
        """
        if self.units == "Jy":
            speed_of_light = 2.99792458e8
            factor = speed_of_light / (self.wavelength**2 * 1e-6 * 1e26)
            self.flux *= factor
            self.flux_err *= factor
            self.units = "W/m2"

    @classmethod
    def from_vizier(
        cls, 
        target: str, 
        radius_arcsec: float = 5.0, 
        aggregate: bool = True
    ) -> "Dataset":
        """
        Fetches SED data from the VizieR SED service.

        Args:
            target (str): Name of the object (resolvable by Simbad/NED).
            radius_arcsec (float): Search radius in arcseconds.
            aggregate (bool): If True, averages multiple points sharing the same filter.

        Returns:
            Dataset: An initialized Dataset object with the fetched data.

        Raises:
            RuntimeError: If the connection to VizieR fails or data is corrupted.
        """
        try:
            url = (f"https://vizier.cds.unistra.fr/viz-bin/sed?"
                   f"-c={target}&-c.rs={radius_arcsec}")
            table = Table.read(url, format='votable')
        except Exception as e:
            raise RuntimeError(f"Failed to fetch data from VizieR for '{target}': {e}")
        
        # Initial cleanup: replace NaN errors with 0.0
        if 'sed_eflux' in table.colnames:
            table['sed_eflux'].fill_value = 0.0
            table['sed_eflux'] = table['sed_eflux'].filled()

        # Aggregate data points by filter if requested
        if aggregate and 'sed_filter' in table.colnames:
            table = cls._aggregate_vizier_table(table)

        # Sort by frequency (standard for SEDs)
        table.sort('sed_freq', reverse=True)

        # Physical conversion: GHz -> um and Jy -> W/m2/um
        freq_ghz = table['sed_freq'].data
        wl_um = 299792.458 / freq_ghz
        
        # Constant: 1e-26 [Jy->W/m2/Hz] * c [m/s] / (lambda [m])^2 / 1e-6 [m->um]
        c = 299792458.0
        factor = 1e-26 * (c / (wl_um * 1e-6)**2) * 1e-6
        flux_wm2 = table['sed_flux'].data * factor
        err_wm2 = table['sed_eflux'].data * factor

        filters = table['sed_filter'].data if 'sed_filter' in table.colnames else None

        return cls(
            wavelength=wl_um, 
            flux=flux_wm2, 
            flux_err=err_wm2, 
            filters=filters, 
            table=table
        )

    @staticmethod
    def _aggregate_vizier_table(table: Table) -> Table:
        """
        Groups points by filter and averages their values.

        Args:
            table (Table): The raw VizieR table.

        Returns:
            Table: A new table with unique filters and averaged flux/frequency.
        """
        grouped = table.group_by('sed_filter')
        new_rows = []
        
        for group in grouped.groups:
            f_mean = np.mean(group['sed_flux'])
            freq_mean = np.mean(group['sed_freq'])
            e_mean = np.mean(group['sed_eflux'])
            
            # Use the first filter name found in the group
            new_rows.append((group['sed_filter'][0], freq_mean, f_mean, e_mean))
        
        return Table(
            rows=new_rows, 
            names=['sed_filter', 'sed_freq', 'sed_flux', 'sed_eflux']
        )
    
    def apply_extinction_correction(self, ebv: float, rv: float = 3.1):
        """
        Applies interstellar de-reddening correction to the dataset.

        Args:
            ebv (float): E(B-V) color excess.
            rv (float): Total-to-selective extinction ratio. Defaults to 3.1.
        """
        self.flux = apply_extinction(self.wavelength, self.flux, ebv, rv)
        if self.flux_err.size > 0:
            self.flux_err = apply_extinction(self.wavelength, self.flux_err, ebv, rv)

    def restrict(self, condition: str):
        """
        Filters data based on a logical condition string.

        Args:
            condition (str): A string like 'wavelength > 1.0' or 'flux < 1e-12'.
                             Evaluated against the internal arrays.
        """
        namespace = {
            'self': self,
            'wavelength': self.wavelength,
            'flux': self.flux,
            'eflux': self.flux_err
        }
        try:
            mask = eval(condition, {"__builtins__": None}, namespace)
            self.wavelength = self.wavelength[mask]
            self.flux = self.flux[mask]
            if self.flux_err.size > 0:
                self.flux_err = self.flux_err[mask]
            if self.filters is not None and self.filters.size > 0:
                self.filters = self.filters[mask]
        except Exception as e:
            print(f"Error applying restriction '{condition}': {e}")

    def plot(self, ax: Optional[plt.Axes] = None, **kwargs) -> plt.Axes:
        """
        Plots the observational SED data.

        Args:
            ax (Optional[plt.Axes]): Existing matplotlib axes. If None, creates a new one.
            **kwargs: Additional keyword arguments passed to ax.errorbar.

        Returns:
            plt.Axes: The matplotlib axes object.
        """
        if ax is None:
            _, ax = plt.subplots()
        
        ax.errorbar(self.wavelength, self.flux, yerr=self.flux_err, fmt='+', **kwargs)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r"Wavelength [$\mu$m]")
        ax.set_ylabel(f"Flux [{self.units}]")
        
        return ax