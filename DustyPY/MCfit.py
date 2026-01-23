import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Any, Optional, List

# Core components from PyMCMC
from PyMCMC import FunctionFitter, Prior, MCMCAnalyzer, ParallelTempering
# Internal DustyPY imports
from .Data import Data  # Correctly import the class to avoid AttributeError
from . import utils as utils

class fit():
    """
    Management class for MCMC fitting procedures in DustyPY.
    
    This class handles the transition between Dusty physical parameters and 
    the PyMCMC inference engine. It supports both single-core and parallel sampling.

    Attributes:
        _Data (Data): The observational data to fit.
        _ParamFit (dict): MCMC configuration (iterations, method, etc.).
        _Param (dict): Physical model parameters and their priors.
        _Results (dict): Dictionary containing the results of the fitting procedure.
        _Stats (dict): Summary statistics for the sampled parameters.
        Ncpu (int): Number of CPUs dedicated to sampling.
        fitter (FunctionFitter): The PyMCMC engine instance.
        analyzer (MCMCAnalyzer): Post-processing tool for diagnostics and plots.
    """

    def __init__(self, data: Optional[Data] = None, ParamFit: Optional[Dict] = None, 
                 Param: Optional[Dict] = None, ncpu: int = 1) -> None:
        """
        Initializes the fitting class.

        Args:
            data (Data, optional): The data object. Defaults to an empty Data instance.
            ParamFit (dict, optional): Sampling settings. Defaults to standard DRAM configuration.
            Param (dict, optional): Model parameters dictionary. Defaults to empty.
            ncpu (int, optional): Number of workers for parallel execution. Defaults to 1.
        """
        # Fix: use Data() directly since it refers to the class imported via __init__.py
        self._Data = data if data is not None else Data()
        
        if ParamFit is None:
            ParamFit = {
                'nsimu': 10000,
                'method': 'DRAM',
                'adaptint': 100,
                'waitbar': True,
                'burn_in_fraction': 0.2
            }
        
        self._ParamFit = ParamFit
        self._Param = Param if Param is not None else {}
        self._Results = {}
        self._Stats = None
        self.Ncpu = ncpu
        self._UserDefinedObject = None
        self.fitter = None   
        self.analyzer = None 

    def set_Data(self, data: Optional[Data] = None, user_defined_object: Optional[list] = None) -> None:
        """Sets the dataset and the contextual object for the physical model."""
        if data is not None:
            self._Data = data
        if user_defined_object is not None:
            self._UserDefinedObject = user_defined_object

    def get_Data(self) -> Data:
        """Returns the current data object."""
        return self._Data

    def set_ParamFit(self, ParamFit: dict) -> None:
        """Updates MCMC simulation settings."""
        self._ParamFit = ParamFit

    def get_ParamFit(self) -> dict:
        """Returns the current MCMC simulation settings."""
        return self._ParamFit

    def set_Param(self, Param: dict) -> None:
        """Updates the physical parameters configuration."""
        self._Param = Param

    def get_Param(self) -> dict:
        """Returns the physical parameters configuration."""
        return self._Param

    def get_Results(self) -> dict:
        """Returns the MCMC chain and results dictionary."""
        return self._Results
    
    def get_Stats(self) -> Any:
        """Returns the summary statistics of the parameters."""
        return self._Stats

    def fit(self, Chi2_func: Any = utils.chi2) -> None:
        """
        Executes the fitting procedure using PyMCMC with a bridge wrapper.
        """
        y = self._Data.get_ydata()
        x = self._Data.get_xdata()
        y_err = self._Data.get_yerr()

        # 1. Map parameters and convert to Priors
        priors, initial_theta0 = utils.convert_params_to_priors(self._Param)
        
        # 2. Context bundle for Dusty
        # Format: [DustyObj, DataObj, FitObj, Logfile, Jansky, Lock]
        data_bundle = [self._UserDefinedObject[0], self._Data, self, False, True, None]

        # 3. BRIDGE WRAPPER : Corrects signature and injects context
        # PyMCMC calls: model(x_input, theta)
        # DustyPY expects: utils.model(theta, context_object)
        def model_bridge(x_input, theta):
            # Create a lightweight context object that mimics what utils.model expects
            class DataContext:
                def __init__(self, x_val, bundle):
                    self.xdata = x_val
                    self.user_defined_object = [bundle]
            
            ctx = DataContext(x_input, data_bundle)
            return utils.model(theta, ctx)

        # 4. LOG-LIKELIHOOD WRAPPER (For the sampler)
        # Note: We use the same context logic for consistency
        def log_likelihood_wrapper(theta):
            class DataContext:
                def __init__(self, x_val, y_val, bundle):
                    self.xdata = [x_val]
                    self.ydata = [y_val]
                    self.user_defined_object = [bundle]
            ctx = DataContext(x, y, data_bundle)
            return -0.5 * Chi2_func(theta, ctx)

        # 5. Initialize the Fitter with the BRIDGE model
        self.fitter = FunctionFitter(
            model_func=model_bridge, # Use the bridge here!
            x_data=x,
            y_data=y,
            y_err=y_err,
            priors=priors,
            custom_log_lik=log_likelihood_wrapper,
            adapt_every=self._ParamFit.get('adaptint', 100)
        )
        self.fitter.user_defined_object = [data_bundle]

        # 6. Execute Sampling
        method = self._ParamFit.get('method', 'DRAM').upper()
        nsimu = self._ParamFit.get('nsimu', 10000)
        show_progress = self._ParamFit.get('waitbar', True)

        if self.Ncpu > 1:
            starts = [np.array(initial_theta0) * (1 + 0.01 * np.random.randn(len(initial_theta0))) 
                      for _ in range(self.Ncpu)]
            chains = self.fitter.fit_parallel(
                initial_params_list=starts,
                n_iterations=nsimu,
                method=method,
                num_workers=self.Ncpu
            )
            self.analyzer = MCMCAnalyzer(chains, fitter=self.fitter, 
                                        burn_in_fraction=self._ParamFit.get('burn_in_fraction', 0.2))
        else:
            chain = self.fitter.fit(
                initial_params=initial_theta0,
                n_iterations=nsimu,
                method=method,
                show_progress=show_progress
            )
            self.analyzer = MCMCAnalyzer(chain, fitter=self.fitter,
                                        burn_in_fraction=self._ParamFit.get('burn_in_fraction', 0.2))

        # 7. Finalize Results
        self._Results = self.analyzer.merged_chain.to_dict()
        self._Stats = self.analyzer.get_summary_stats()

    def print_results(self, include_ic: bool = True) -> None:
        """
        Prints the summary statistics.
        Set include_ic=False if the model is slow to skip AIC/BIC/WAIC.
        """
        if self.analyzer:
            param_names = [k for k, v in self._Param.items() if v['sample']]
            self.analyzer.print_summary(param_names=param_names, include_ic=include_ic)
            
    def plot_stats(self) -> None:
        """Plots parameter traces and posterior distributions (Corner Plot)."""
        if self.analyzer:
            self.analyzer.plot_traces()
            self.analyzer.plot_corner()

    def plot_pairwise_correlation(self, fig: Optional[dict] = None) -> None:
        """Generates a Corner Plot to visualize parameter correlations."""
        if self.analyzer:
            self.analyzer.plot_corner()

    def prediction_interval(self, x_eval: np.ndarray, n_samples: int = 500) -> Dict[str, np.ndarray]:
        """
        Calculates confidence intervals for the model across a custom grid.

        Args:
            x_eval (np.ndarray): The wavelength grid for evaluation.
            n_samples (int): Number of posterior samples to draw.
        
        Returns:
            dict: Statistical bounds (mean, median, 68%, 95%).
        """
        samples = self.analyzer.merged_chain.samples
        indices = np.random.choice(len(samples), size=min(n_samples, len(samples)), replace=False)
        
        y_models = []
        for idx in indices:
            theta = samples[idx]
            # Use the prediction model utility
            y_models.append(utils.prediction_model(theta, self.fitter.user_defined_object[0], x_eval))
            
        y_models = np.array(y_models)
        
        return {
            'mean': np.mean(y_models, axis=0),
            'median': np.median(y_models, axis=0),
            'low_95': np.percentile(y_models, 2.5, axis=0),
            'high_95': np.percentile(y_models, 97.5, axis=0),
            'low_68': np.percentile(y_models, 16, axis=0),
            'high_68': np.percentile(y_models, 84, axis=0)
        }

    def plot_prediction_interval(self, wavelength_grid: np.ndarray, fig_settings: Optional[Dict] = None) -> None:
        """Plots the SED with prediction intervals overlaid on observational data."""
        intervals = self.prediction_interval(wavelength_grid)
        
        plt.figure(figsize=(10, 6))
        plt.fill_between(wavelength_grid, intervals['low_95'], intervals['high_95'], color='gray', alpha=0.3, label='95% CI')
        plt.fill_between(wavelength_grid, intervals['low_68'], intervals['high_68'], color='gray', alpha=0.5, label='68% CI')
        plt.plot(wavelength_grid, intervals['median'], color='red', label='Median Model')
        
        plt.errorbar(self._Data.get_xdata(), self._Data.get_ydata(), yerr=self._Data.get_yerr(), fmt='ok', label='Data')
        
        if fig_settings:
            plt.xscale(fig_settings.get('xscale', 'log'))
            plt.yscale(fig_settings.get('yscale', 'log'))
            plt.xlabel(fig_settings.get('xlabel', 'Wavelength ($\mu m$)'))
            plt.ylabel(fig_settings.get('ylabel', 'Flux'))
            
        plt.legend()
        plt.show()