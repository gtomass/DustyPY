try:
    from . import utils as utils
    from .constants import *
except ImportError:
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
    import utils as utils
    from constants import *
import numpy as np

class Data():
    """
    Class representing data with x and y values and their respective errors.
    """

    def __init__(self, xdata: np.array = None, ydata: np.array = None, xerr: np.array = None, yerr: np.array = None) -> None:
        """
        Initializes an instance of the Data class.

        Parameters:
        xdata (array-like, optional): The x data values. Defaults to None.
        ydata (array-like, optional): The y data values. Defaults to None.
        xerr (array-like, optional): The errors in the x data values. Defaults to None.
        yerr (array-like, optional): The errors in the y data values. Defaults to None.
        """
        self._xdata = xdata
        self._ydata = ydata
        self._xerr = xerr
        self._yerr = yerr

    def set_xdata(self, xdata: np.array = None) -> None:
        """
        Sets the x data values. 

        Parameters:
        xdata (array-like): The x data values.
        """
        self._xdata = xdata

    def set_ydata(self, ydata: np.array = None) -> None:
        """
        Sets the y data values.

        Parameters:
        ydata (array-like): The y data values.
        """
        self._ydata = ydata

    def set_xerr(self, xerr: np.array = None) -> None:
        """
        Sets the errors in the x data values.

        Parameters:
        xerr (array-like): The errors in the x data values.
        """
        self._xerr = xerr

    def set_yerr(self, yerr: np.array = None) -> None:
        """
        Sets the errors in the y data values.

        Parameters:
        yerr (array-like): The errors in the y data values.
        """
        self._yerr = yerr

    def get_xdata(self) -> np.array:
        """
        Returns the x data values.

        Returns:
        array-like: The x data values.
        """
        return self._xdata
    
    def get_ydata(self) -> np.array:
        """
        Returns the y data values.

        Returns:
        array-like: The y data values.
        """
        return self._ydata
    
    def get_xerr(self) -> np.array:
        """
        Returns the errors in the x data values.

        Returns:
        array-like: The errors in the x data values.
        """
        return self._xerr

    def get_yerr(self) -> np.array:
        """
        Returns the errors in the y data values.

        Returns:
        array-like: The errors in the y data values.
        """
        return self._yerr
    
    def ImportData(self, Path: str) -> np.array:
        """
        Imports data from a specified file path.

        Parameters:
        Path (str): The file path to import data from.

        Returns:
        array-like: The imported data.
        """
        extention = Path.split('.')[-1]

        if extention == 'csv':
            return utils.LoadCSV(Path)
        elif extention in ['fits', 'fit']:
            return utils.LoadFits(Path)
        else:
            return utils.LoadFile(Path)

    def ScatterData(self, 
                    unit: str = None, 
                    xlim: tuple = None, 
                    ylim: tuple = None, 
                    ax = None, 
                    scale: str = 'linear', 
                    kwargs: dict = None) -> None:
        """
        Plots the data as a scatter plot or error plot.

        Parameters:
        unit (str, optional): The unit of the axes. Defaults to None.
        xlim (tuple, optional): The limits of the x-axis. Defaults to None.
        ylim (tuple, optional): The limits of the y-axis. Defaults to None.
        ax (matplotlib.axes.Axes, optional): The axis on which to plot. Defaults to None.
        scale (str, optional): The scale of the axes ('linear' or 'log'). Defaults to 'linear'.
        kwargs (dict, optional): Additional arguments for the plot function. Defaults to None.
        """
        if self._yerr is None:
            utils.ScatterPlot(self._ydata, self._xdata, unit=unit, xlim=xlim, ylim=ylim, ax=ax, scale=scale, kwargs=kwargs)
        else:
            utils.ErrorPlot(self._ydata, self._xdata, self._yerr, unit=unit, xlim=xlim, ylim=ylim, ax=ax, scale=scale, kwargs=kwargs)

    def UnredData(self, EBV: float, Rv:float = 3.1) -> None:
        """
        Applies a dereddening correction to the data.

        Parameters:
        EBV (float): The E(B-V) value for the correction.
        Rv (float, optional): The ratio of total to selective extinction. Defaults to 3.1.

        Warning:
        Wavelength must be in micrometers. 
        """
        self._ydata = utils.Unred(self._xdata * 1e4, self._ydata, EBV=EBV, Rv=Rv)
        self._yerr = utils.Unred(self._xdata * 1e4, self._yerr, EBV=EBV, Rv=Rv) if self._yerr is not None else None

    def ConvertJansky(self) -> None:
        """
        Converts the y data and its errors from Watts per square meter to Jansky.
        """
        self._ydata = utils.WattToJansky(self._ydata, self._xdata * um)
        self._yerr = utils.WattToJansky(self._yerr, self._xdata * um) if self._yerr is not None else None
            
    def ConvertWatt(self) -> None:
        """
        Converts the y data and its errors from Jansky to Watts per square meter.
        """
        self._ydata = utils.JanskyToWatt(self._ydata, self._xdata * um)
        self._yerr = utils.JanskyToWatt(self._yerr, self._xdata * um) if self._yerr is not None else None

    def QuerryVizierData(self, radius: float = 5, target: str = None):
        """
        Queries Vizier for data.

        Parameters:
        radius (str): The radius of the query.
        target (str): The target of the query.

        Returns:
        array-like: The queried data.
        """
        if target is None:
            raise ValueError('target must be specified')
        return utils.VizierQuery(radius, target)
    
    def SetVizierData(self, table) -> None:
        """
        Sets the data from a Vizier query. Fill nan with 0 for compatibility with the rest of the code.

        Parameters:
        table (array-like): The data from the Vizier query.
        """
        self._xdata = np.nan_to_num(np.asarray(c / (table['sed_freq'] * Ghz) * 1/um),nan=0.)
        self._ydata = np.nan_to_num(np.asarray(table['sed_flux']),nan=0.)
        self._yerr = np.nan_to_num(np.asarray(table['sed_eflux']),nan=0.)
    
    def RestrictData(self, ListOfCondition = list[str]):
        """
        Restricts the data based on a condition.
        Parameters:
        condition (str): The condition to restrict the data with.

        Examples:
        >>> data = Data(xdata=np.array([1, 2, 3, 4]), ydata=np.array([10, 20, 30, 40]))
        >>> data.RestrictData(['xdata > 2'])
        >>> print(data.get_xdata())
        [3 4]
        >>> print(data.get_ydata())
        [30 40]
        """
        for condition in ListOfCondition:
            condition = 'self._'+condition.split(' ')[0] + ' ' + ' '.join(condition.split(' ')[1:])
            restriction = eval(condition)

            self._xdata = self._xdata[restriction]
            self._ydata = self._ydata[restriction]
            self._xerr = self._xerr[restriction] if self._xerr is not None else None
            self._yerr = self._yerr[restriction] if self._yerr is not None else None