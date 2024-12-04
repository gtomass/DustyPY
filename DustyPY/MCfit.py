import pymcmcstat
import pymcmcstat.MCMC
from . import utils as utils
from . import Data as Data
import numpy as np


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

    def __init__(self, data: Data.Data = None, Model: pymcmcstat.MCMC.MCMC = None, ParamFit: dict = None, Param: dict = None) -> None:
        """
        Initializes an instance of the Fit class.

        Parameters:
        Data (Data, optional): The data to fit. Defaults to an instance of Data.
        Model (pymcmcstat.MCMC.MCMC, optional): The MCMC model. Defaults to an instance of pymcmcstat.MCMC.MCMC.
        ParamFit (dict, optional): The fitting parameters. Defaults to a predefined dictionary.
        Param (dict, optional): The model parameters. Defaults to an empty dictionary.
        """
        if data is None:
            data = Data.Data()
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
        self._Stats = None

    def set_Data(self, data: Data.Data = None) -> None:
        """
        Sets the data to fit.

        Parameters:
        Data (Data, optional): The data to fit. Defaults to None.
        """
        if data is None:
            data = Data.Data()
        self._Data = data

    def get_Data(self) -> Data.Data:
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

    def fit(self, Chi2=utils.chi2) -> None:
        """
        Performs the fitting procedure using the provided chi-squared function.

        Parameters:
        Chi2 (function, optional): The chi-squared function to be used for the fitting procedure. Defaults to utils.chi2.
        """
        y = self._Data.get_ydata()
        x = self._Data.get_xdata()
        self._Model.data.add_data_set(x, y)
        self.set_Model()

        self.set_Chi2Func(Chi2)
        self._Model.run_simulation()

        results = {}
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


    # def PredictionModel(self) -> None:
    #     """
    #     Sets up the prediction model for calculating prediction intervals.
    #     """
    #     nds = len(self._Data.get_xdata())
    #     print(nds)

    #     def pred_modelfun(preddata, theta):
    #         return utils.model(theta[:2], preddata.xdata[0], preddata.ydata[0]).reshape(nds,)

    #     self._Model.PI.setup_prediction_interval_calculation(pred_modelfun)
