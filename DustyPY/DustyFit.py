from . import dusty
from . import MCfit
from . import Data
from . import utils as utils
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import time
import pymcmcstat.MCMC


class DustyFit():
    """
    Class representing a fitting procedure for a Dusty model.

    Attributes:
    _Dusty (Dusty): The Dusty model used to fit the data.
    _Data (Data): The data to fit.
    _ParamFit (dict): The fitting parameters.
    _Param (dict): The model parameters.
    _Fit (fit): The fitting object.
    """

    def __init__(self, Dusty: dusty.Dusty, data: Data.Data, ParamFit: dict = None, Param: dict = None, fit: MCfit.fit = None, logfile: bool = False) -> None:
        """
        Initializes an instance of the DustyFit class.

        Parameters:
        Dusty (Dusty): The Dusty model to be fitted.
        data (Data): The data to fit.
        ParamFit (dict, optional): The fitting parameters. Defaults to None.
        Param (dict, optional): The model parameters. Defaults to None.
        fit (fit, optional): The fitting object. Defaults to an instance of MCfit.fit.
        """
        if fit is None:
            fit = MCfit.fit()
        if data is None:
            data = Data.Data()
        

        self._Dusty = Dusty
        self._Data = data
        self._ParamFit = ParamFit
        self._Param = Param
        self._Fit = fit
        self._logfile = logfile

    def __InitFit(self, Jansky: bool = True) -> None:
        """
        Initializes the fitting object with the data, fitting parameters, and model parameters.
        """
        self._Fit.set_Data(data=self._Data, user_defined_object = [self._Dusty, self._Data, self._Fit, self._logfile, Jansky])
        self._Fit.set_ParamFit(ParamFit=self._ParamFit)
        self._Fit.set_Param(Param=self._Param)

    # Méthodes set
    def set_Dusty(self, Dusty: dusty.Dusty) -> None:
        """
        Sets the Dusty model to be fitted.

        Parameters:
        Dusty (Dusty): The Dusty model to be fitted.
        """
        self._Dusty = Dusty
    

    def set_Data(self, Data: Data.Data) -> None:
        """
        Sets the data to fit.

        Parameters:
        Data (Data): The data to fit.
        """
        self._Data = Data

    def set_ParamFit(self, ParamFit: dict) -> None:
        """
        Sets the fitting parameters.

        Parameters:
        ParamFit (dict): The fitting parameters.
        """
        self._ParamFit = ParamFit

    def set_Param(self, Param: dict) -> None:
        """
        Sets the model parameters.

        Parameters:
        Param (dict): The model parameters.
        """
        self._Param = Param

    def set_Fit(self, fit: MCfit.fit) -> None:
        self._Fit = fit

    # Méthodes get
    def get_Dusty(self) -> dusty.Dusty:
        """
        Returns the Dusty model to be fitted.

        Returns:
        Dusty: The Dusty model to be fitted.
        """
        return self._Dusty

    def get_Data(self) -> Data.Data:
        """
        Returns the data to fit.

        Returns:
        Data: The data to fit.
        """
        return self._Data

    def get_ParamFit(self) -> dict:
        """
        Returns the fitting parameters.

        Returns:
        dict: The fitting parameters.
        """
        return self._ParamFit

    def get_Param(self) -> dict:
        """
        Returns the model parameters.

        Returns:
        dict: The model parameters.
        """
        return self._Param

    def get_Fit(self) -> MCfit.fit:
        """
        Returns the fitting object.

        Returns:
        fit: The fitting object.
        """
        return self._Fit

    def __setChange(self, change: dict) -> None:
        """
        Applies changes to the Dusty model based on the provided dictionary.

        Parameters:
        change (dict): A dictionary containing the changes to be applied. The keys should specify the parameter to change (e.g., 'Temp', 'Lum', 'Opacity' for each stars of the model) and the values should be the new values for those parameters.

        Raises:
        NotImplementedError: If an attempt is made to change the dust size, which is not yet fittable.
        """
        utils.set_change(dusty=self._Dusty,change=change)
            

    def __Chi2Dusty(self, theta, data) -> float:
        """
        Calculate the chi-squared value for the Dusty model fit.
        This method updates the Dusty model parameters based on the provided 
        theta values, runs the Dusty model, and computes the chi-squared value 
        between the model's spectral energy distribution (SED) and the observed data.
        Parameters:
        -----------
        theta : list or array-like
            A list or array of parameter values to update the Dusty model.
        data : object
            An object containing the observed data with attributes `xdata` and `ydata`.
            If `xdata` and `ydata` are not attributes, the object should have methods 
            `get_xdata()` and `get_ydata()` to retrieve the data.
        Returns:
        --------
        float
            The chi-squared value representing the goodness of fit between the model 
            SED and the observed data. If observational errors (`yerr`) are available, 
            they are used in the chi-squared calculation.
        """

        try:
            ydata = data.ydata[0]
        except AttributeError:
            ydata = data.get_ydata()

        ymodel = utils.model(theta, data).reshape(ydata.shape)

        chi2 = np.nansum(((ymodel - ydata)/self._Data.get_yerr())**2) if self._Data.get_yerr() is not None else np.nansum((ymodel - ydata)**2)

        return chi2
        
    def __Chi2Dusty_modified(self, theta, data) -> float:
        """
        Calculate the chi-squared value for the Dusty model fit.
        This method updates the Dusty model parameters based on the provided 
        theta values, runs the Dusty model, and computes the chi-squared value 
        between the model's spectral energy distribution (SED) and the observed data.
        Parameters:
        -----------
        theta : list or array-like
            A list or array of parameter values to update the Dusty model.
        data : object
            An object containing the observed data with attributes `xdata` and `ydata`.
            If `xdata` and `ydata` are not attributes, the object should have methods 
            `get_xdata()` and `get_ydata()` to retrieve the data.
        Returns:
        --------
        float
            The chi-squared value representing the goodness of fit between the model 
            SED and the observed data. If observational errors (`yerr`) are available, 
            they are used in the chi-squared calculation.
        """

        try:
            xdata = data.xdata[0]
            ydata = data.ydata[0]
        except AttributeError:
            xdata, ydata = data.get_xdata(), data.get_ydata()

        fdata_ks = ydata[np.argmin(abs(xdata-2.190))]
        ymodel = utils.model(theta, data).reshape(ydata.shape)

        if np.isnan(ymodel).all():
            return np.inf

        fmodel_ks = ymodel[np.argmin(abs(xdata-2.190))]
        ymodel_norm = ymodel/fmodel_ks
        ydata_norm = ydata/fdata_ks

        chi2 = 1/(len(ydata)-len(theta)-1)*np.nansum((1-(ymodel_norm/ydata_norm))**2/(ymodel_norm/ydata_norm))


        return chi2

    def lunch_fit(self, chi2: str = 'Chi2', logfile: bool = False, Jansky: bool = True) -> None:
        """
        Initializes the fitting procedure and performs the fit using the chi-squared function.

        This method initializes the fitting object and then runs the fitting procedure using the chi-squared function specific to the Dusty model.
        """
        self._logfile = logfile
        begin = time.time()

        self.__InitFit(Jansky=Jansky)
        if chi2 == 'Chi2':
            self._Fit.fit(Chi2=self.__Chi2Dusty)
        elif chi2 == 'Chi2_modified':
            self._Fit.fit(Chi2=self.__Chi2Dusty_modified)
        else:
            raise Exception('The chi2 function is not recognized')
        
        print(f'Fitting time: {time.time()-begin} s')

    def get_chi2(self, chi2: str = 'Chi2', Jansky: bool = True) -> float:
        """
        Returns the chi-squared value of the fitting procedure.

        Returns: chi2 (float): The chi-squared value of the fitting procedure.
        """
        pdata = pymcmcstat.MCMC.MCMC()
        pdata.data.add_data_set(self.get_Data().get_xdata(), self.get_Data().get_ydata(), user_defined_object=[self._Dusty, self._Data, self._Fit, self._logfile, Jansky])

        p = [key for key, value in self.get_Param().items() if value['sample']]
        stats = [self._Fit.get_Stats()['mean'][p.index(key)] if value['sample'] else value['theta0'] for i, (key, value) in enumerate(self.get_Param().items())]
        if chi2 == 'Chi2':
            return self.__Chi2Dusty(stats, pdata.data)
        elif chi2 == 'Chi2_modified':
            return self.__Chi2Dusty_modified(stats, pdata.data)
        else:
            raise Exception('The chi2 function is not recognized')

    def print_results(self) -> None:
        """
        Prints the results of the fitting procedure.

        This method calls the print_results method of the fitting object to display the results of the fitting procedure.
        """
        self._Fit.print_results()

    def plot_stats(self) -> None:
        """
        Plots the statistics of the fitting procedure.

        Parameters:
        ax (matplotlib.axes.Axes): The axis on which to plot the statistics.
        """
        self._Fit.plot_stats()

    def plot_results(self,
                     unit: dict = None,
                     xlim: tuple = None,
                     ylim: tuple = None,
                     ax: plt.Axes = None,
                     scale: str = 'linear',
                     Jansky: bool = True,
                     kwargs_fit: dict = None,
                     kwargs_data: dict = None,
                     normalize: bool = False,
                     SED_band: bool = False,
                     save: bool = False) -> None:
        """
        Plots the results of the fitting procedure along with the data.

        Parameters:
        unit (str, optional): What is plotted on each axes. Defaults to None.
        xlim (tuple, optional): The limits of the x-axis. Defaults to None.
        ylim (tuple, optional): The limits of the y-axis. Defaults to None.
        ax (matplotlib.axes.Axes, optional): The axis on which to plot. Defaults to None.
        scale (str, optional): The scale of the axes ('linear' or 'log'). Defaults to 'linear'.
        kwargs_fit (dict, optional): Additional arguments for the fit plot. Defaults to None.
        kwargs_data (dict, optional): Additional arguments for the data plot. Defaults to None.

        This method retrieves the fitting results, updates the Dusty model parameters, runs the Dusty simulation, and generates the Spectral Energy Distribution (SED). It then plots the SED and the data with optional error bars.
        """
        sampled_params = [key for key, value in self._Fit.get_Param().items() if value['sample']]
        try:
            Lum = self._Fit.get_Stats()['mean'][sampled_params.index('Lest')]
        except:
            Lum = self._Fit.get_Param()['Lest']['theta0']

        stats_mean = self._Fit.get_Stats()['mean']
        dustsize = self._Dusty.get_Model().get_Dust().get_DustSize()
        density = self._Dusty.get_Model().get_Dust().get_Density()

        if 'amin' in sampled_params:
            dustsize['amin'] = np.round(10**stats_mean[sampled_params.index('amin')],3)
        if 'amax' in sampled_params:
            dustsize['amax'] = np.round(10**stats_mean[sampled_params.index('amax')],3)
        if 'q' in sampled_params:
            dustsize['q'] = stats_mean[sampled_params.index('q')]
        if 'shell' in sampled_params:
            density['shell'] = stats_mean[sampled_params.index('shell')]

        results_params = dict(zip(sampled_params, stats_mean))
        for param in ['amin', 'amax', 'q', 'shell']:
            if param in results_params.keys():
                results_params.pop(param, None)        
                
        results_params.update({'DustSize':dustsize})
        results_params.update({'Density':density})
        self.__setChange(change=results_params)
        self._Dusty.change_parameter()
        self._Dusty.lunch_dusty(verbose = 0)
        self._Dusty.make_SED(distance=self._Dusty.get_Model().get_Distance(), luminosity=Lum, Jansky=Jansky)
        SED = self._Dusty.get_SED()
        utils.plot(SED.get_Flux(), SED.get_Wavelength(), unit=unit,
                   xlim=xlim, ylim=ylim, ax=ax, scale=scale, kwargs=kwargs_fit, normalize=normalize)
        
        if SED_band:
            SED.scatter_SED_bandpass(self._Data.get_common_filters(self._Data.get_table()), ax=ax, kwargs={
                                     'marker': '.', 'color': 'b', 'label': 'SED'}, normalize=normalize)
        if self._Data.get_yerr() is not None:
            utils.error_plot(self._Data.get_ydata(), self._Data.get_xdata(), self._Data.get_yerr(
            ), unit=unit, xlim=xlim, ylim=ylim, ax=ax, scale=scale, kwargs=kwargs_data, normalize=normalize)
        else:
            if 'fmt' in kwargs_data.keys():
                marker = kwargs_data['fmt']
                kwargs_data.pop('fmt')
                kwargs_data.update({'marker': marker})
            utils.scatter_plot(self._Data.get_ydata(), self._Data.get_xdata(
            ), unit=unit, xlim=xlim, ylim=ylim, ax=ax, scale=scale, kwargs=kwargs_data, normalize=normalize)

        if save:
            plt.savefig(os.path.join(self._Dusty.get_PATH(), self._Dusty.get_Model().get_Name(),'SED.png'), dpi=300)

        plt.legend()
        plt.show()


    def plot_interval(self, wavelength_intervals: list = [(1e-2, 1e0, 50), (1e0, 1e2, 100), (1e2, 1e4, 50)],ciset: dict = None, piset: dict = None, fig: tuple = None) -> None:
        """
        Plots the interval of the fitting procedure.
        """
        wavelengths = []
        for start, end, num in wavelength_intervals:
            wavelengths.extend(utils.log_space(start, end, num))
        
        wavelengths = np.asarray(wavelengths)
        self._Fit.plot_prediction_interval(wavelength_dusty = wavelengths, ciset = ciset, piset = piset, fig = fig)

    def plot_pairwise_correlation(self, fig: dict = None) -> None:
        """
        Plots the pairwise correlation of the fitting procedure.
        """
        self._Fit.plot_pairwise_correlation(fig = fig)