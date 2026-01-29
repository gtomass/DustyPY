"""
Likelihood and objective function module for DustyPY.
Defines various statistical methods to compare model outputs with observational data.
"""

from enum import Enum, auto
from abc import ABC, abstractmethod
from typing import Optional, Union, Callable, Any

import numpy as np


class LikelihoodType(Enum):
    """Enumeration of supported likelihood/chi-square types."""
    STANDARD = auto()
    SHAPE = auto()
    CUSTOM = auto()


class Likelihood(ABC):
    """
    Abstract Base Class for likelihood calculations.
    
    Any objective function used by the Fitter must inherit from this class 
    and implement the 'compute' method.
    """

    @abstractmethod
    def compute(
        self, 
        y_obs: np.ndarray, 
        y_model: np.ndarray, 
        y_err: np.ndarray, 
        n_params: int
    ) -> float:
        """
        Calculates the cost/likelihood value.

        Args:
            y_obs (np.ndarray): Observed flux values.
            y_model (np.ndarray): Model-predicted flux values.
            y_err (np.ndarray): Observational uncertainties.
            n_params (int): Number of free parameters in the fit.

        Returns:
            float: The calculated cost value (e.g., reduced chi-square).
        """
        pass


class StandardChi2(Likelihood):
    """
    Implements a standard reduced chi-square calculation.
    
    Compares absolute flux values. If uncertainties are zero or negative, 
    a default 10% error floor is applied to prevent division by zero.
    """

    def compute(
        self, 
        y_obs: np.ndarray, 
        y_model: np.ndarray, 
        y_err: np.ndarray, 
        n_params: int
    ) -> float:
        """
        Calculates the reduced chi-square for absolute flux values.

        Args:
            y_obs: Observed fluxes.
            y_model: Model fluxes.
            y_err: Observed flux errors.
            n_params: Number of free parameters.

        Returns:
            float: Reduced chi-square value. Returns np.inf if degrees of 
                freedom are non-positive.
        """
        # Degrees of Freedom (DOF)
        dof = len(y_obs) - n_params - 1
        if dof <= 0:
            return np.inf

        # Safety: Apply 10% error floor for points with no defined error
        err = np.where(y_err > 0, y_err, 0.1 * y_obs)
        
        chi2 = np.nansum(((y_model - y_obs) / err)**2)
        return chi2 / dof


class ShapeChi2(Likelihood):
    """
    Implements a normalized shape-based chi-square calculation.
    
    This method normalizes both the model and observations at a specific 
    wavelength (pivot) before comparison. It is useful when fitting the 
    spectral energy distribution (SED) shape independently of absolute 
    luminosity or distance.
    """

    def __init__(self, wavelengths: np.ndarray, norm_wavelength: float = 2.19):
        """
        Initializes the ShapeChi2 objective function.

        Args:
            wavelengths (np.ndarray): Wavelength grid associated with the fluxes.
            norm_wavelength (float): Pivot wavelength used for normalization (um).
                Defaults to 2.19 um (K-band).
        """
        self.wavelengths = wavelengths
        self.norm_wavelength = norm_wavelength

    def compute(
        self, 
        y_obs: np.ndarray, 
        y_model: np.ndarray, 
        y_err: np.ndarray, 
        n_params: int
    ) -> float:
        """
        Calculates the cost based on normalized SED shape.

        Args:
            y_obs: Observed fluxes.
            y_model: Model fluxes.
            y_err: Observed flux errors (ignored in current shape logic).
            n_params: Number of free parameters.

        Returns:
            float: Normalized shape cost. Returns np.inf if normalization 
                fails or if degrees of freedom are invalid.
        """
        # Find index closest to the normalization wavelength
        idx_norm = np.argmin(np.abs(self.wavelengths - self.norm_wavelength))
        
        # Check if normalization is valid
        if y_model[idx_norm] <= 0 or y_obs[idx_norm] <= 0:
            return np.inf
        
        # Normalize both model and observations
        y_mod_norm = y_model / y_model[idx_norm]
        y_obs_norm = y_obs / y_obs[idx_norm]
        
        dof = len(y_obs) - n_params - 1
        if dof <= 0:
            return np.inf
        
        # Calculate ratio-based deviation
        ratio = y_mod_norm / y_obs_norm
        # Standardize result by degrees of freedom
        return (1.0 / dof) * np.nansum((1.0 - ratio)**2 / ratio)


def get_likelihood(
    l_type: Union[str, LikelihoodType], 
    dataset: Optional[Any] = None, 
    custom_func: Optional[Likelihood] = None, 
    **kwargs
) -> Union[Likelihood, Callable]:
    """
    Factory function to retrieve the appropriate likelihood instance.

    Args:
        l_type (Union[str, LikelihoodType]): The type of likelihood to retrieve.
            Accepts strings like 'standard' or 'shape'.
        dataset (Optional[Dataset]): Required for 'SHAPE' type to access 
            wavelengths.
        custom_func (Optional[Likelihood]): An instance of a custom Likelihood 
            class.
        **kwargs: Additional parameters passed to the Likelihood constructor 
            (e.g., norm_wavelength).

    Returns:
        Likelihood: An instance of the requested likelihood class.

    Raises:
        ValueError: If an unknown type is requested or if required 
            arguments are missing.
    """
    # Convert string to Enum member if necessary
    if isinstance(l_type, str):
        try:
            l_type = LikelihoodType[l_type.upper()]
        except KeyError:
            raise ValueError(f"Unknown likelihood string identifier: {l_type}")

    if l_type == LikelihoodType.STANDARD:
        return StandardChi2()
    
    elif l_type == LikelihoodType.SHAPE:
        if dataset is None:
            raise ValueError("The 'SHAPE' likelihood requires a dataset to "
                             "access wavelength information.")
        return ShapeChi2(wavelengths=dataset.wavelength, **kwargs)
    
    elif l_type == LikelihoodType.CUSTOM:
        if custom_func is None:
            raise ValueError("A custom likelihood object must be provided "
                             "when using LikelihoodType.CUSTOM.")
        return custom_func
    
    raise ValueError(f"Unsupported likelihood type: {l_type}")