import pymcmcstat
import pymcmcstat.MCMC
import pymcmcstat.propagation
from . import utils as utils
from . import Data as Data
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool, Lock

def run_simulation_wrapper(Model):
        """
        Wrapper function to run the MCMC simulation.
        """
        Model.run_simulation()


class fit():
    """
    Class representing a fitting procedure using MCMC.

    Attributes:
    _Model (pymcmcstat.MCMC.MCMC): The MCMC model.
    _Data (Data): The data to fit.
    _ParamFit (dict): The fitting parameters.
    _Param (dict): The model parameters.
    _Results (dict): The results of the fitting procedure.
    """

    def __init__(self, data: Data = None, Model: pymcmcstat.MCMC.MCMC = None, ParamFit: dict = None, Param: dict = None, ncpu: int = 1) -> None:
        """
        Initializes an instance of the Fit class.

        Parameters:
        Data (Data, optional): The data to fit. Defaults to an instance of Data.
        Model (pymcmcstat.MCMC.MCMC, optional): The MCMC model. Defaults to an instance of pymcmcstat.MCMC.MCMC.
        ParamFit (dict, optional): The fitting parameters. Defaults to a predefined dictionary.
        Param (dict, optional): The model parameters. Defaults to an empty dictionary.
        """
        if data is None:
            data = Data()
        if Model is None:
            Model = pymcmcstat.MCMC.MCMC()
        if ParamFit is None:
            ParamFit = {
                'nsimu': 10000,
                'updatesigma': True,
                'method': 'dram',
                'adaptint': 100,
                'verbosity': 0,
                'waitbar': True,
            }
        if Param is None:
            Param = {}

        self._Model = Model
        self._Data = data
        self._ParamFit = ParamFit
        self._Param = Param
        self._Results = {}
        self.Ncpu = ncpu
        self._UserDefinedObject = None
        self._Stats = None

    def set_Data(self, data: Data = None, user_defined_object: list = None) -> None:
        """
        Sets the data to fit.

        Parameters:
        Data (Data, optional): The data to fit. Defaults to None.
        """
        if data is None:
            data = Data()
        if user_defined_object is not None:
            self._UserDefinedObject = user_defined_object
        self._Data = data

    def get_Data(self) -> Data:
        """
        Returns the data to fit.

        Returns:
        Data: The data to fit.
        """
        return self._Data

    def set_Model(self) -> None:
        """
        Sets the MCMC model parameters and simulation options.

        Utilizes the utility function `set_mcmc_param` to set the model parameters and
        defines the simulation options using the fitting parameters.
        """
        utils.set_mcmc_param(self._Model, self._Param)
        self._Model.simulation_options.define_simulation_options(
            **self._ParamFit)

    def get_Model(self) -> pymcmcstat.MCMC.MCMC:
        """
        Returns the MCMC model.

        Returns:
        pymcmcstat.MCMC.MCMC: The MCMC model.
        """
        return self._Model

    def set_ParamFit(self, ParamFit: dict = None) -> None:
        """
        Sets the fitting parameters for the MCMC model.

        Parameters:
        ParamFit (dict, optional): A dictionary containing the fitting parameters. If None, defaults to a predefined dictionary:
        {
            'nsimu': 10000,
            'updatesigma': True,
            'method': 'dram',
            'adaptint': 100,
            'verbosity': 0,
            'waitbar': True,
        }
        """
        if ParamFit is None:
            raise ValueError('ParamFit cannot be None')
        self._ParamFit = ParamFit

    def get_ParamFit(self) -> dict:
        """
        Returns the fitting parameters for the MCMC model.

        Returns:
        dict: The fitting parameters.
        """
        return self._ParamFit

    def set_Param(self, Param: dict) -> None:
        """
        Sets the model parameters for the MCMC model.

        Parameters:
        Param (dict): A dictionary containing the model parameters.
        """
        if Param is None:
            raise ValueError('Param cannot be None')
        self._Param = Param

    def get_Param(self) -> dict:
        """
        Returns the model parameters for the MCMC model.

        Returns:
        dict: The model parameters.
        """
        return self._Param

    def set_Chi2Func(self, func=None) -> None:
        """
        Sets the chi-squared function for the MCMC model.

        Parameters:
        func (function): The chi-squared function to be used by the MCMC model.
        """
        if func is None:
            raise ValueError('func cannot be None')
        self._Model.model_settings.define_model_settings(sos_function=func)

    def get_Results(self) -> dict:
        """
        Returns the results of the fitting procedure.

        Returns:
        dict: The results of the fitting procedure.
        """
        return self._Results
    
    def get_Stats(self) -> any:
        """
        Returns the statistics of the fitting procedure.

        Returns:
        any: The statistics of the fitting procedure.
        """
        return self._Stats
    
    def run_parallel(self) -> None:
        """
        Runs the MCMC simulation in parallel using multiple processes.
        """
        self._Model.run_simulation()


    def fit(self, Chi2=utils.chi2) -> None:
        """
        Executes the fitting procedure using the specified chi-squared function.

        Parameters:
        Chi2 (function, optional): The chi-squared function to use for the fitting procedure. Defaults to utils.chi2.
        """
        y = self._Data.get_ydata()
        x = self._Data.get_xdata()

        self.set_Chi2Func(Chi2)

        # Add data and user-defined object to the model
        user_defined_object = self._UserDefinedObject
        if self.Ncpu > 1:
            lock = Lock()
            user_defined_object.append(lock)
        else:
            user_defined_object.append(None)

        self._UserDefinedObject = user_defined_object
        self._Model.data.add_data_set(x, y, user_defined_object=self._UserDefinedObject)
        self.set_Model()

        if self.Ncpu > 1:
            with Pool(processes=self.Ncpu) as pool:
                from functools import partial
                pool.map(partial(run_simulation_wrapper, self._Model), range(self.Ncpu))
        else:
            self._Model.run_simulation()

        # Store results and statistics
        results = self._Model.simulation_results.results.copy()
        chain = results['chain']
        burnin = int(results['nsimu'] / 2)
        stats = self._Model.chainstats(chain[burnin:, :], results=results, returnstats=True)

        self._Results = results
        self._Stats = stats

    def print_results(self) -> None:
        """
        Prints the results of the fitting procedure, including chain statistics.
        """
        chain = self._Results['chain']
        burnin = int(self._Results['nsimu'] / 2)
        self._Model.chainstats(chain[burnin:, :], self._Results)

    def plot_stats(self) -> None:
        """
        Plots the statistics of the fitting procedure.

        Parameters:
        ax (matplotlib.axes.Axes): The axes to plot the statistics on.
        """
        chain = self.get_Results()['chain']
        burnin = int(self.get_Results()['nsimu'] / 2)
        names = self.get_Results()['names']
        s2chain = self.get_Results()['s2chain']
        
        result = self.get_Model().chainstats(
            chain[burnin:, :], self.get_Results(), returnstats=True)
        mcpl = pymcmcstat.mcmcplot # initialize plotting methods
        mcpl.plot_density_panel(chain[burnin:,:], names);
        mcpl.plot_chain_panel(chain[burnin:,:], names);
        mcpl.plot_density_panel(np.sqrt(s2chain[burnin:,:]), ['$\\sigma_1$', '$\\sigma_2$'])

    def plot_density_panel(self) -> None:
        """
        Plots the density panel of the MCMC chain.
        """
        chain = self.get_Results()['chain']
        burnin = int(self.get_Results()['nsimu'] / 2)
        names = self.get_Results()['names']
        mcpl = pymcmcstat.mcmcplot
        mcpl.plot_density_panel(chain[burnin:,:], names)

    def plot_chain_panel(self) -> None:
        """
        Plots the density panel of the MCMC chain.
        """
        chain = self.get_Results()['chain']
        burnin = int(self.get_Results()['nsimu'] / 2)
        names = self.get_Results()['names']
        mcpl = pymcmcstat.mcmcplot
        mcpl.plot_chain_panel(chain[burnin:,:], names)


    def prediction_interval(self, data: any = None) -> None:
        """
        Sets up the prediction model for calculating prediction intervals.
        """
        def predmodelfun2(theta, dat):
            p = utils.prediction_model(theta,dat)
            return p

        results = self.get_Results()
        chain = results['chain']
        s2chain = results['s2chain']

        intervals = pymcmcstat.propagation.calculate_intervals(chain, results, data, predmodelfun2,
                               s2chain=s2chain, nsample=500, waitbar=True)
        
        return intervals
    
    def plot_prediction_interval(self, wavelength_dusty: np.array = None, ciset: dict = None, piset: dict = None, fig: dict = None) -> None:
        """
        Plots the prediction intervals.

        Parameters:
        intervals (dict): The prediction intervals to plot.
        """
        def format_plot(fig):
            plt.xscale(fig['xscale'])
            plt.yscale(fig['yscale'])
            plt.xlim(*fig['xlim'])
            plt.ylim(*fig['ylim'])
            plt.xlabel(fig['xlabel'], fontsize=15)
            plt.ylabel(fig['ylabel'], fontsize=15)
            plt.title(fig['title'])

        if wavelength_dusty is None:
            time = self._Model.data.xdata[0]
        else:
            pdata = pymcmcstat.MCMC.MCMC()
            pdata.data.add_data_set(wavelength_dusty, wavelength_dusty, user_defined_object=self._UserDefinedObject)
            time = wavelength_dusty
        if ciset is None:
            ciset = {'limits': [50, 95]}
        if piset is None:
            addprediction = False
        else:
            addprediction = True

        intervals = self.prediction_interval(data=pdata.data)


        if fig is None:
            pymcmcstat.propagation.plot_intervals(intervals=intervals, time = time, 
                                              ydata=self._Data.get_ydata(), xdata=self._Data.get_xdata(), adddata=True,
                                              ciset=ciset, piset=piset, addprediction=addprediction)
        else:
            f,ax = pymcmcstat.propagation.plot_intervals(intervals=intervals, time = time, 
                                                    ydata=self._Data.get_ydata(), xdata=self._Data.get_xdata(), adddata=True,
                                                    ciset=ciset, piset=piset, addprediction=addprediction)
            format_plot(fig)
        

    def plot_pairwise_correlation(self, fig: dict = None) -> None:
        """
        Plots the pairwise correlation between the parameters.
        """
        chain = self.get_Results()['chain']
        names = self.get_Results()['names']
        settings = dict(fig=fig)
        mcpl = pymcmcstat.mcmcplot.plot_pairwise_correlation_panel(chain, names, settings=settings)