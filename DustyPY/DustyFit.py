from . import dusty 
from . import MCfit
from . import Data
from . import utils as utils
import numpy as np
import matplotlib.pyplot as plt
import subprocess

class DustyFit():
    """
    Class representing a fitting procedure for a Dusty model.

    Attributes:
    _Dusty (Dusty): The Dusty model used to fit the data.
    _Data (Data): The data to fit.
    _ParamFit (dict): The fitting parameters.
    _Param (dict): The model parameters.
    _Fit (Fit): The fitting object.
    """

    def __init__(self, Dusty, Data, ParamFit=None, Param=None, Fit=MCfit.Fit()):
        """
        Initializes an instance of the DustyFit class.

        Parameters:
        Dusty (Dusty): The Dusty model to be fitted.
        Data (Data): The data to fit.
        ParamFit (dict, optional): The fitting parameters. Defaults to None.
        Param (dict, optional): The model parameters. Defaults to None.
        Fit (Fit, optional): The fitting object. Defaults to an instance of MCfit.Fit.
        """
        self._Dusty = Dusty
        self._Data = Data
        self._ParamFit = ParamFit
        self._Param = Param
        self._Fit = Fit

    def __InitFit(self):
        """
        Initializes the fitting object with the data, fitting parameters, and model parameters.
        """
        self._Fit.set_Data(Data=self._Data)
        self._Fit.set_ParamFit(ParamFit=self._ParamFit)
        self._Fit.set_Param(Param=self._Param)

    # Méthodes set
    def set_Dusty(self, Dusty):
        """
        Sets the Dusty model to be fitted.

        Parameters:
        Dusty (Dusty): The Dusty model to be fitted.
        """
        self._Dusty = Dusty

    def set_Data(self, Data):
        """
        Sets the data to fit.

        Parameters:
        Data (Data): The data to fit.
        """
        self._Data = Data

    def set_ParamFit(self, ParamFit):
        """
        Sets the fitting parameters.

        Parameters:
        ParamFit (dict): The fitting parameters.
        """
        self._ParamFit = ParamFit

    def set_Param(self, Param):
        """
        Sets the model parameters.

        Parameters:
        Param (dict): The model parameters.
        """
        self._Param = Param

    def set_Fit(self, Fit):
        self._Fit = Fit

    # Méthodes get
    def get_Dusty(self):
        """
        Returns the Dusty model to be fitted.

        Returns:
        Dusty: The Dusty model to be fitted.
        """
        return self._Dusty

    def get_Data(self):
        """
        Returns the data to fit.

        Returns:
        Data: The data to fit.
        """
        return self._Data

    def get_ParamFit(self):
        """
        Returns the fitting parameters.

        Returns:
        dict: The fitting parameters.
        """
        return self._ParamFit

    def get_Param(self):
        """
        Returns the model parameters.

        Returns:
        dict: The model parameters.
        """
        return self._Param

    def get_Fit(self):
        """
        Returns the fitting object.

        Returns:
        Fit: The fitting object.
        """
        return self._Fit
    
    def __setChange(self, change):
        """
        Applies changes to the Dusty model based on the provided dictionary.

        Parameters:
        change (dict): A dictionary containing the changes to be applied. The keys should specify the parameter to change (e.g., 'Temp', 'Lum', 'Opacity' for each stars of the model) and the values should be the new values for those parameters.

        Raises:
        NotImplementedError: If an attempt is made to change the dust size, which is not yet fittable.
        """
        for key in change.keys():
            if 'Temp' in key:
                self._Dusty.get_Model().get_Stars()[int(key.split('Temp')[-1])-1].set_Temperature(change[key])
            elif 'Lum' in key:
                self._Dusty.get_Model().get_Stars()[int(key.split('Lum')[-1])-1].set_Luminosity(change[key])
            elif key in ['Opacity']:
                self._Dusty.get_Model().get_Dust().set_tau(change[key])
            elif key in ['Composition', 'Abundances']:
                self._Dusty.get_Model().get_Dust().set_Composition(change[key])
            elif key in ['DustSize']:
                self._Dusty.get_Model().get_Dust().set_DustSize(change[key])
            elif key == 'Lest':
                pass
            else:
                raise NotImplementedError(f'Parameter {key} not recognized.')
            
    def __Chi2Dusty(self,theta,data):
        
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
        
        change = utils.ListToDict(list(self._Fit.get_Param().keys()),theta)
        

        self.__setChange(change)
        
        self._Dusty.ChangeParameter()
        self._Dusty.LunchDusty(verbose = 0)
        self._Dusty.MakeSED(distance = self._Dusty.get_Model().get_Distance())

        try:
            xdata = data.xdata[0]
            ydata = data.ydata[0]
        except AttributeError:
            xdata,ydata = data.get_xdata(),data.get_ydata()

        ymodel = utils.model(theta[-1], xdata, self._Dusty.GetSED().get_Wavelength(), self._Dusty.GetSED().get_Flux()
                             ).reshape(ydata.shape)
        
        subprocess.call('clear', shell=True)

        if self._Data.get_yerr() is not None:
            return np.nansum(((ymodel - ydata)/self._Data.get_yerr())**2)
        else:    
            return np.nansum((ymodel - ydata)**2)
            
    
    def LunchFit(self):
        """
        Initializes the fitting procedure and performs the fit using the chi-squared function.

        This method initializes the fitting object and then runs the fitting procedure using the chi-squared function specific to the Dusty model.
        """
        self.__InitFit()
        self._Fit.Fit(Chi2=self.__Chi2Dusty)

    def PrintResults(self):
        """
        Prints the results of the fitting procedure.

        This method calls the PrintResults method of the fitting object to display the results of the fitting procedure.
        """
        self._Fit.PrintResults()
    
    def PlotResults(self, unit=None, xlim=None, ylim=None, ax=None, scale='linear', kwargs_fit=None, kwargs_data=None):
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

        chain = self._Fit.get_Results()['chain']
        burnin = int(self._Fit.get_Results()['nsimu'] / 2)
        result = self._Fit.get_Model().chainstats(chain[burnin:, :], self._Fit.get_Results(), returnstats=True)

        self.__setChange(utils.ListToDict(list(self._Fit.get_Param().keys()), result['mean']))
        self._Dusty.ChangeParameter()
        self._Dusty.LunchDusty()
        self._Dusty.MakeSED(distance=self._Dusty.get_Model().get_Distance())
        SED = self._Dusty.GetSED()
        utils.Plot(SED.get_Flux(), SED.get_Wavelength(), unit=unit, xlim=xlim, ylim=ylim, ax=ax, scale=scale, kwargs=kwargs_fit)
        if self._Data.get_yerr() is not None:
            utils.ErrorPlot(self._Data.get_ydata(), self._Data.get_xdata(), self._Data.get_yerr(), unit=unit, xlim=xlim, ylim=ylim, ax=ax, scale=scale, kwargs=kwargs_data)
        else:
            if 'fmt' in kwargs_data.keys():
                marker = kwargs_data['fmt']
                kwargs_data.pop('fmt')
                kwargs_data.update({'marker': marker})
            utils.ScatterPlot(self._Data.get_ydata(), self._Data.get_xdata(), unit=unit, xlim=xlim, ylim=ylim, ax=ax, scale=scale, kwargs=kwargs_data)
        plt.show()


    def __Check(self):
        pass
