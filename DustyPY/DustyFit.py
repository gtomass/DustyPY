import os
import time
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Optional, Any, List, Union

from . import dusty
from . import MCfit
from . import Data
from . import utils as utils

class DustyFit():
    """
    High-level class to manage the fitting of a Dusty model to observational data.

    This class acts as a bridge between the physical Dusty model, the observational 
    dataset, and the MCMC sampling engine provided by PyMCMC.

    Attributes:
        _Dusty (dusty.Dusty): The Dusty model instance used to generate synthetic spectra.
        _Data (Data.Data): The observational data to be fitted.
        _ParamFit (dict): Settings for the MCMC simulation.
        _Param (dict): Definition of physical parameters (initial values, bounds, etc.).
        _Fit (MCfit.fit): The sampling manager utilizing the PyMCMC library.
        _logfile (bool): Whether to enable Dusty's internal logging.
    """

    def __init__(self, Dusty: dusty.Dusty, data: Data.Data, ParamFit: Optional[dict] = None, 
                 Param: Optional[dict] = None, fit_obj: Optional[MCfit.fit] = None, 
                 logfile: bool = False, ncpu: int = 1) -> None:
        """
        Initializes the DustyFit procedure.

        Args:
            Dusty (dusty.Dusty): The Dusty model instance.
            data (Data.Data): The observational data.
            ParamFit (dict, optional): MCMC simulation settings.
            Param (dict, optional): Model parameter configurations.
            fit_obj (MCfit.fit, optional): Pre-existing fit object. Defaults to a new MCfit.fit instance.
            logfile (bool, optional): Enable Dusty logging. Defaults to False.
            ncpu (int, optional): Number of CPU cores for parallel sampling. Defaults to 1.
        """
        if fit_obj is None:
            fit_obj = MCfit.fit(ncpu=ncpu)
        if data is None:
            data = Data.Data()
        
        self._Dusty = Dusty
        self._Data = data
        self._ParamFit = ParamFit
        self._Param = Param
        self._Fit = fit_obj
        self._logfile = logfile

    def __InitFit(self, Jansky: bool = True) -> None:
        """
        Prepares the internal fit object with context and parameters before execution.
        """
        # User defined object bundle passed to the sampler
        context = [self._Dusty, self._Data, self._Fit, self._logfile, Jansky]
        self._Fit.set_Data(data=self._Data, user_defined_object=context)
        self._Fit.set_ParamFit(ParamFit=self._ParamFit)
        self._Fit.set_Param(Param=self._Param)

    # Setter methods
    def set_Dusty(self, Dusty: dusty.Dusty) -> None:
        """Sets the Dusty model instance."""
        self._Dusty = Dusty

    def set_Data(self, Data: Data.Data) -> None:
        """Sets the observational data."""
        self._Data = Data

    def set_ParamFit(self, ParamFit: dict) -> None:
        """Sets the MCMC simulation parameters."""
        self._ParamFit = ParamFit

    def set_Param(self, Param: dict) -> None:
        """Sets the physical model parameters."""
        self._Param = Param

    def set_Fit(self, fit_obj: MCfit.fit) -> None:
        """Sets the sampling manager object."""
        self._Fit = fit_obj

    # Getter methods
    def get_Dusty(self) -> dusty.Dusty:
        """Returns the Dusty model instance."""
        return self._Dusty

    def get_Data(self) -> Data.Data:
        """Returns the observational data."""
        return self._Data

    def get_ParamFit(self) -> dict:
        """Returns the MCMC settings."""
        return self._ParamFit

    def get_Param(self) -> dict:
        """Returns the model parameters."""
        return self._Param

    def get_Fit(self) -> MCfit.fit:
        """Returns the sampling manager."""
        return self._Fit

    def __setChange(self, change: dict) -> None:
        """
        Applies parameter changes to the underlying Dusty model.
        """
        utils.set_change(dusty=self._Dusty, change=change)

    def Chi2Dusty(self, theta: list, data: Any) -> float:
        """
        Calculates the standard Chi-squared value for the model fit.
        
        Args:
            theta (list): Current parameter vector being sampled.
            data: Data context provided by the sampler.
            
        Returns:
            float: Reduced chi-squared value.
        """
        try:
            ydata = data.ydata[0]
        except AttributeError:
            ydata = data.get_ydata()

        try:
            ymodel = utils.model(theta, data).reshape(ydata.shape)
        except ValueError as e:
            raise ValueError(f"Model computation failed: {e}")

        # Compute chi2 using observational errors if available
        if self._Data.get_yerr() is not None:
            chi2 = np.nansum(((ymodel - ydata) / self._Data.get_yerr())**2)
        else:
            chi2 = np.nansum((ymodel - ydata)**2)

        # Degrees of freedom calculation
        dof = len(ydata) - len(theta) - 1
        return chi2 / dof if dof > 0 else np.inf

    def Chi2Dusty_modified(self, theta: list, data: Any) -> float:
        """
        Calculates a modified Chi-squared value normalized at the Ks band (2.19 um).
        
        This is often used to ensure the spectral shape is fitted correctly regardless 
        of absolute scaling.
        """
        try:
            xdata = data.xdata[0]
            ydata = data.ydata[0]
        except AttributeError:
            xdata, ydata = data.get_xdata(), data.get_ydata()

        # Find Ks band flux for normalization
        fdata_ks = ydata[np.argmin(abs(xdata - 2.190))]
        ymodel = utils.model(theta, data).reshape(ydata.shape)

        if np.isnan(ymodel).all():
            return np.inf

        fmodel_ks = ymodel[np.argmin(abs(xdata - 2.190))]
        ymodel_norm = ymodel / fmodel_ks
        ydata_norm = ydata / fdata_ks

        dof = len(ydata) - len(theta) - 1
        chi2 = (1 / dof) * np.nansum((1 - (ymodel_norm / ydata_norm))**2 / (ymodel_norm / ydata_norm))
        return chi2

    def lunch_fit(self, chi2: str = 'Chi2', logfile: bool = False, Jansky: bool = True) -> None:
        """
        Launches the MCMC fitting procedure.

        Args:
            chi2 (str): Identifier for the cost function ('Chi2' or 'Chi2_modified').
            logfile (bool): Enable Dusty logging.
            Jansky (bool): Use Jansky units for flux.
        """
        self._logfile = logfile
        start_time = time.time()

        self.__InitFit(Jansky=Jansky)
        
        if chi2 == 'Chi2':
            self._Fit.fit(Chi2_func=self.Chi2Dusty)
        elif chi2 == 'Chi2_modified':
            self._Fit.fit(Chi2_func=self.Chi2Dusty_modified)
        else:
            raise ValueError(f"Cost function '{chi2}' is not recognized.")
        
        print(f"Fitting completed in: {time.time() - start_time:.2f} s")

    def get_chi2(self, chi2: str = 'Chi2', Jansky: bool = True) -> float:
        """
        Calculates the Chi-squared value for the current parameter means.
        """
        # Create a lightweight context object for the chi2 calculation
        # This replaces the legacy pymcmcstat requirement
        class Context:
            def __init__(self, x, y, bundle):
                self.xdata = [x]
                self.ydata = [y]
                self.user_defined_object = [bundle]

        bundle = [self._Dusty, self._Data, self._Fit, self._logfile, Jansky, None]
        ctx = Context(self.get_Data().get_xdata(), self.get_Data().get_ydata(), bundle)

        # Identify sampled parameters and retrieve their mean values from the Analyzer
        p_keys = [key for key, value in self.get_Param().items() if value['sample']]
        stats = []
        for key, value in self.get_Param().items():
            if value['sample']:
                idx = p_keys.index(key)
                stats.append(self._Fit.get_Stats()['mean'][idx])
            else:
                stats.append(value['theta0'])

        if chi2 == 'Chi2':
            return self.Chi2Dusty(stats, ctx)
        elif chi2 == 'Chi2_modified':
            return self.Chi2Dusty_modified(stats, ctx)
        else:
            raise ValueError(f"Cost function '{chi2}' not recognized.")

    def print_results(self) -> None:
        """Prints the statistical results of the fit."""
        self._Fit.print_results()

    def plot_stats(self) -> None:
        """Plots chain traces and parameter distributions."""
        self._Fit.plot_stats()

    def plot_results(self, unit: dict = None, xlim: tuple = None, ylim: tuple = None, 
                     ax: plt.Axes = None, scale: str = 'linear', Jansky: bool = True, 
                     kwargs_fit: dict = None, kwargs_data: dict = None, 
                     normalize: bool = False, SED_band: bool = False, 
                     save: bool = False) -> None:
        """
        Plots the fitted SED against the observational data.
        """
        sampled_params = [key for key, value in self._Fit.get_Param().items() if value['sample']]
        stats_mean = self._Fit.get_Stats()['mean']

        # Extract Luminosity (Lest)
        try:
            Lum = stats_mean[sampled_params.index('Lest')]
        except (ValueError, IndexError):
            Lum = self._Param.get('Lest', {}).get('theta0', 1e4)

        # Extract and update dust properties
        dustsize = self._Dusty.get_Model().get_Dust().get_DustSize()
        density = self._Dusty.get_Model().get_Dust().get_Density()

        if 'amin' in sampled_params:
            dustsize['amin'] = np.round(10**stats_mean[sampled_params.index('amin')], 3)
        if 'amax' in sampled_params:
            dustsize['amax'] = np.round(10**stats_mean[sampled_params.index('amax')], 3)
        if 'q' in sampled_params:
            dustsize['q'] = stats_mean[sampled_params.index('q')]
        if 'shell' in sampled_params:
            density['shell'] = stats_mean[sampled_params.index('shell')]

        # Update and run Dusty simulation
        results_params = dict(zip(sampled_params, stats_mean))
        for p in ['amin', 'amax', 'q', 'shell']:
            results_params.pop(p, None)        
                
        results_params.update({'DustSize': dustsize, 'Density': density})
        self.__setChange(change=results_params)
        self._Dusty.change_parameter()
        self._Dusty.lunch_dusty(verbose=0)
        self._Dusty.make_SED(luminosity=Lum, Jansky=Jansky)
        
        # Plotting
        SED = self._Dusty.get_SED()
        utils.plot(SED.get_Flux(), SED.get_Wavelength(), unit=unit, xlim=xlim, 
                   ylim=ylim, ax=ax, scale=scale, kwargs=kwargs_fit, normalize=normalize)
        
        if SED_band:
            SED.scatter_SED_bandpass(self._Data.get_common_filters(), ax=ax, 
                                     kwargs={'fmt': '.b', 'label': 'SED Filters'}, 
                                     normalize=normalize)
        
        # Plot data points
        if self._Data.get_yerr() is not None:
            utils.error_plot(self._Data.get_ydata(), self._Data.get_xdata(), 
                             self._Data.get_yerr(), unit=unit, xlim=xlim, ylim=ylim, 
                             ax=ax, scale=scale, kwargs=kwargs_data, normalize=normalize)
        else:
            utils.scatter_plot(self._Data.get_ydata(), self._Data.get_xdata(), unit=unit, 
                               xlim=xlim, ylim=ylim, ax=ax, scale=scale, kwargs=kwargs_data, 
                               normalize=normalize)

        if save:
            save_path = os.path.join(self._Dusty.get_PATH(), self._Dusty.get_Model().get_Name(), 'SED.png')
            plt.savefig(save_path, dpi=300)

        plt.legend()

    def plot_interval(self, wavelength_intervals: list = [(1e-2, 1e0, 50), (1e0, 1e2, 100), (1e2, 1e4, 50)],
                      ciset: dict = None, piset: dict = None, fig_settings: dict = None) -> None:
        """
        Evaluates and plots the posterior prediction intervals across a wavelength grid.
        """
        wavelengths = []
        for start, end, num in wavelength_intervals:
            wavelengths.extend(utils.log_space(start, end, num))
        
        self._Fit.plot_prediction_interval(wavelength_grid=np.asarray(wavelengths), 
                                          fig_settings=fig_settings)

    def plot_pairwise_correlation(self, fig: dict = None) -> None:
        """Plots the parameter correlation matrix (Corner Plot)."""
        self._Fit.plot_pairwise_correlation(fig=fig)