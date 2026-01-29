"""
MCMC Fitting module for DustyPY.
Bridges DustyPY with the PyMCMC library to perform parameter estimation 
using Markov Chain Monte Carlo methods.
"""

import os
import shutil
from typing import Dict, List, Optional, Any, Union

import numpy as np
from PyMCMC import FunctionFitter, Prior, MCMCAnalyzer

from .likelihood import LikelihoodType, get_likelihood


class Fitter:
    """
    Orchestrates the MCMC fitting process for Dusty models.

    This class manages the mapping between MCMC sampler variables and 
    physical Dusty parameters. It handles parallel process isolation, 
    physical constraint checking, and likelihood computation.

    Attributes:
        model (Model): The DustyPY model instance to fit.
        runner (Runner): The execution engine for Dusty.
        dataset (Dataset): Observational data to fit against.
        n_workers (int): Number of CPU cores for parallel sampling.
        params_config (Dict): Configuration of free parameters.
        likelihood (Likelihood): Likelihood function object.
    """

    def __init__(self, model, runner, dataset, n_workers: int = 1):
        """
        Initializes the Fitter.

        Args:
            model (Model): The model object to optimize.
            runner (Runner): The runner instance to execute Dusty calls.
            dataset (Dataset): The observational data.
            n_workers (int): Number of parallel workers. Defaults to 1.
        """
        self.model = model
        self.runner = runner
        self.dataset = dataset
        self.n_workers = n_workers

        self.params_config: Dict[str, Any] = {}
        self.likelihood = None
        self._current_fit_type = "photometry"

    def _setup_prior(self, prior_type: str, **kwargs) -> Prior:
        """
        Configures a PyMCMC Prior object based on provided arguments.

        Args:
            prior_type (str): Type of prior ('uniform', 'gaussian', etc.).
            **kwargs: Parameters for the prior (e.g., min_val, max_val, mu, sigma).

        Returns:
            Prior: A configured Prior instance.

        Raises:
            ValueError: If the prior type is unsupported or arguments are missing.
        """
        if prior_type == 'uniform':
            params = [kwargs.get('min_val'), kwargs.get('max_val')]
        elif prior_type in ['gaussian', 'normal']:
            params = [kwargs.get('mu'), kwargs.get('sigma')]
        elif prior_type == 'log_normal':
            params = [kwargs.get('mu'), kwargs.get('sigma')]
        elif prior_type == 'beta':
            params = [kwargs.get('alpha'), kwargs.get('beta')]
        else:
            raise ValueError(f"Prior type '{prior_type}' is not supported.")

        if any(p is None for p in params):
            raise ValueError(
                f"Missing required arguments for prior '{prior_type}'. "
                f"Received: {list(kwargs.keys())}"
            )

        return Prior(prior_type, params)

    def add_parameter(self, name: str, initial: float, prior: str = 'uniform', **kwargs):
        """
        Adds a global system parameter to the fit.

        Args:
            name (str): Parameter name (e.g., 'opacity', 'distance', 'rgd').
            initial (float): Starting value for the sampler.
            prior (str): Prior distribution type.
            **kwargs: Arguments for the prior distribution.
        """
        prior_obj = self._setup_prior(prior, **kwargs)
        self.params_config[name] = {
            'initial': initial,
            'prior': prior_obj,
            'type': 'global'
        }

    def add_star_parameter(self, star_index: int, star_prop: str, initial: float,
                           prior: str = 'uniform', **kwargs):
        """
        Adds a stellar parameter to the fit.

        Args:
            star_index (int): Index of the star in the model's star list.
            star_prop (str): Property to fit ('temp', 'lum', or 'logg').
            initial (float): Starting value.
            prior (str): Prior distribution type.
            **kwargs: Arguments for the prior distribution.

        Raises:
            ValueError: If star_prop is not valid.
        """
        if star_prop not in ['temp', 'lum', 'logg']:
            raise ValueError("star_prop must be 'temp', 'lum', or 'logg'")

        param_id = f"star{star_index}_{star_prop}"
        self.params_config[param_id] = {
            'type': 'star',
            'star_index': star_index,
            'star_prop': star_prop,
            'initial': initial,
            'prior': self._setup_prior(prior, **kwargs)
        }

    def _update_model(self, theta: np.ndarray):
        """
        Injects values from the sampler into the Model object.

        Args:
            theta (np.ndarray): Array of parameter values proposed by the sampler.
        """
        for i, (name, cfg) in enumerate(self.params_config.items()):
            val = theta[i]

            if cfg['type'] == 'star':
                idx = cfg['star_index']
                # Ensure the star list is initialized
                while len(self.model.stars) <= idx:
                    self.model.add_star(temp=3000, lum=1.0, logg=4.0)
                prop = cfg['star_prop']
                self.model.stars[idx][prop] = val

            else:  # Global parameters
                if name == 'opacity':
                    self.model.dust['opacity'] = val
                elif name == 'temperature':
                    self.model.dust['temperature'] = val
                elif name == 'power':
                    self.model.geometry['power'] = val
                elif name == 'thickness':
                    self.model.geometry['thickness'] = val
                elif name in ['amin', 'amax', 'q']:
                    sd = self.model.dust.get('size_distribution')
                    if sd:
                        sd[name] = val
                elif name == 'distance':
                    self.model.distance = val
                elif name == 'rgd':
                    self.model.gas_to_dust_ratio = val
                elif name == 'density':
                    self.model.dust_density = val
                else:
                    raise ValueError(f"Unknown global parameter name: {name}")

    def _check_physical_constraints(self, theta: np.ndarray) -> bool:
        """
        Ensures sampler values are within physical and grid-related bounds.

        Args:
            theta (np.ndarray): Proposed parameter vector.

        Returns:
            bool: True if parameters are valid, False otherwise.
        """
        for i, (name, cfg) in enumerate(self.params_config.items()):
            val = theta[i]

            if cfg['type'] == 'star':
                # Temperatures and luminosities must be strictly positive
                if val <= 0:
                    return False
                if cfg['star_prop'] == 'temp' and val < 100.0:
                    return False
                
                # Check atmosphere grid bounds if applicable
                if (cfg['star_prop'] in ['temp', 'logg'] and 
                        self.model.spectral_shape == "atmosphere"):
                    from ..utils.physics import get_grid_bounds
                    t_bounds, g_bounds = get_grid_bounds(self.model.grid_path)
                    if t_bounds and cfg['star_prop'] == 'temp':
                        if not (t_bounds[0] <= val <= t_bounds[1]):
                            return False
                    if g_bounds and cfg['star_prop'] == 'logg':
                        if not (g_bounds[0] <= val <= g_bounds[1]):
                            return False

            else:
                # Most physical parameters must be positive
                positive_params = [
                    'opacity', 'temperature', 'distance', 
                    'rgd', 'amin', 'amax', 'density', 'thickness'
                ]
                if name in positive_params and val <= 0:
                    return False

                if name == 'temperature' and val < 10.0:
                    return False

        return True

    def log_likelihood(self, theta: np.ndarray) -> float:
        """
        Computes the log-likelihood for a given parameter set.

        This method isolates each process execution in a temporary 
        directory to prevent file collision during parallel runs.

        Args:
            theta (np.ndarray): The parameter values.

        Returns:
            float: The log-likelihood value. Returns -inf for invalid points.
        """
        if not self._check_physical_constraints(theta):
            return -np.inf

        fit_type = getattr(self, "_current_fit_type", "photometry")
        self._update_model(theta)

        # Process isolation using PID and random ID
        run_id = f"mcmc_p{os.getpid()}_{np.random.randint(0, 1000000)}"
        tmp_dir = os.path.abspath(os.path.join("tmp_mcmc", run_id))

        try:
            # Runner executes Dusty and computes synthetic photometry
            result = self.runner.run(self.model, run_dir=tmp_dir)

            if fit_type == "photometry":
                y_model, y_obs, y_err = [], [], []
                # Map synthetic photometry to observational data points
                for viz_name, syn_data in result.photometry.items():
                    # Find closest data point in wavelength
                    idx = np.argmin(np.abs(self.dataset.wavelength - syn_data['wavelength']))
                    y_model.append(syn_data['flux'])
                    y_obs.append(self.dataset.flux[idx])
                    
                    # Error handling: use 10% floor if error is missing
                    err = self.dataset.flux_err[idx]
                    if err <= 0:
                        err = 0.1 * self.dataset.flux[idx]
                    y_err.append(err)
                
                y_model, y_obs, y_err = map(np.array, [y_model, y_obs, y_err])
            else:
                # Interpolate the full model spectrum
                y_model = np.interp(self.dataset.wavelength, result.wavelength, result.flux)
                y_obs = self.dataset.flux
                y_err = np.where(self.dataset.flux_err > 0, 
                                 self.dataset.flux_err, 
                                 0.1 * self.dataset.flux)

            # Likelihood computation (e.g., -0.5 * chi2)
            cost = self.likelihood.compute(y_obs, y_model, y_err, len(theta))
            return -0.5 * cost

        except Exception:
            return -np.inf
        finally:
            # Cleanup temporary execution directory
            if os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir, ignore_errors=True)

    def run(
        self,
        n_iterations: int = 1000,
        likelihood: str = 'standard',
        custom_chi2: Optional[Any] = None,
        fit_type: str = "photometry",
        method: str = 'DRAM',
        progress: bool = True,
        threads_per_worker: int = 1,
        **kwargs
    ) -> MCMCAnalyzer:
        """
        Executes the MCMC fit.

        Args:
            n_iterations (int): Number of sampler iterations.
            likelihood (str): Likelihood function name.
            custom_chi2: Optional custom chi2 function.
            fit_type (str): 'photometry' or 'spectrum'.
            method (str): MCMC sampling algorithm (e.g., 'DRAM', 'MH').
            progress (bool): Whether to show the progress bar.
            threads_per_worker (int): Threads per sampler worker.
            **kwargs: Extra arguments for likelihood or sampler.

        Returns:
            MCMCAnalyzer: Analysis object containing chains and statistics.
        """
        self._current_fit_type = fit_type

        initial_theta = [cfg['initial'] for cfg in self.params_config.values()]
        priors = [cfg['prior'] for cfg in self.params_config.values()]

        self.likelihood = get_likelihood(
            l_type=likelihood,
            dataset=self.dataset,
            custom_func=custom_chi2,
            **kwargs
        )

        # Prevent log(0) issues in samplers
        y_err_safe = np.where(self.dataset.flux_err > 0, self.dataset.flux_err, 1e-30)

        fitter = FunctionFitter(
            model_func=None,
            x_data=self.dataset.wavelength,
            y_data=self.dataset.flux,
            y_err=y_err_safe,
            priors=priors,
            custom_log_lik=self.log_likelihood
        )

        if self.n_workers > 1:
            # Generate spread initial positions for parallel chains
            starts = [
                np.array(initial_theta) * (1 + 0.01 * np.random.randn(len(initial_theta)))
                for _ in range(self.n_workers)
            ]
            chains = fitter.fit_parallel(
                starts,
                n_iterations=n_iterations,
                method=method,
                num_workers=self.n_workers,
                show_progress=progress,
                threads_per_worker=threads_per_worker
            )
            analyzer = MCMCAnalyzer(chains, fitter=fitter)
        else:
            chain = fitter.fit(
                initial_theta,
                n_iterations=n_iterations,
                method=method,
                show_progress=progress,
                threads_per_worker=threads_per_worker
            )
            analyzer = MCMCAnalyzer(chain, fitter=fitter)

        # Update the model with the best-fit (mean) values and cache results
        stats = analyzer.get_summary_stats()
        self._update_model(stats['mean'])
        self.best_model = self.runner.run(self.model)

        return analyzer