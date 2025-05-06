import sys
try:
    from . import utils as utils
    from .constants import *
except ImportError:
    import os
    sys.path.append(os.path.abspath(os.path.join(
        os.path.dirname(__file__), os.path.pardir)))
    import utils as utils
    from constants import *
import numpy as np
import astropy.units as u
from astropy.table import unique
from io import StringIO


class Data():
    """
    Class representing data with x and y values and their respective errors.
    """

    def __init__(self, xdata: np.array = None, ydata: np.array = None, xerr: np.array = None, yerr: np.array = None, table = None) -> None:
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
        self._table = table

    def set_data(self, xdata: np.array = None, ydata: np.array = None, xerr: np.array = None, yerr: np.array = None, table = None) -> None:
        """
        Sets the data values.

        Parameters:
        xdata (array-like): The x data values.
        ydata (array-like): The y data values.
        xerr (array-like, optional): The errors in the x data values. Defaults to None.
        yerr (array-like, optional): The errors in the y data values. Defaults to None.
        """
        self._xdata = xdata
        self._ydata = ydata
        self._xerr = xerr
        self._yerr = yerr
        self._table = table

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
    
    def get_table(self):
        """
        Retrieve the table data.
        Returns:
            object: The table data stored in the instance.
        """

        return self._table

    def import_data(self, Path: str, header: int = 0, delimiter: str = ' ') -> np.array:
        """
        Imports data from a specified file path.

        Parameters:
        Path (str): The file path to import data from.
        header (int, optional): The number of header lines to skip. Defaults to 0.
        delimiter (str, optional): The delimiter of the data. Defaults to ' '.

        Returns:
        array-like: The imported data.
        """
        extention = Path.split('.')[-1]

        if extention == 'csv':
            return utils.load_csv(Path, sep=delimiter)
        elif extention in ['fits', 'fit']:
            return utils.load_fits(Path)
        else:
            return self.str_to_data(utils.load_file(Path, header=header))

    def str_to_data(self, data: list) -> np.array:
        """
        Converts a string to data.

        Returns:
        array-like: The data.
        """
        return utils.str_to_data(data)

    def scatter_data(self,
                    unit: str = None,
                    xlim: tuple = None,
                    ylim: tuple = None,
                    ax=None,
                    scale: str = 'linear',
                    kwargs: dict = None,
                    normalize: bool = False) -> None:
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
            marker = kwargs['fmt'] if 'fmt' in kwargs else '+'
            kwargs.pop('fmt', None)
            kwargs.update({'marker': marker})
            
            utils.scatter_plot(self._ydata, self._xdata, unit=unit,
                               xlim=xlim, ylim=ylim, ax=ax, scale=scale, kwargs=kwargs, normalize=normalize)
        else:
            utils.error_plot(self._ydata, self._xdata, self._yerr, unit=unit,
                             xlim=xlim, ylim=ylim, ax=ax, scale=scale, kwargs=kwargs, normalize=normalize)

    def unred_data(self, EBV: float, Rv: float = 3.1) -> None:
        """
        Applies a dereddening correction to the data.

        Parameters:
        EBV (float): The E(B-V) value for the correction.
        Rv (float, optional): The ratio of total to selective extinction. Defaults to 3.1.

        Warning:
        Wavelength must be in micrometers. 
        """
        self._ydata = utils.unred(
            self._xdata * 1e4, self._ydata, EBV=EBV, Rv=Rv)
        self._yerr = utils.unred(
            self._xdata * 1e4, self._yerr, EBV=EBV, Rv=Rv) if self._yerr is not None else None

    def convert_to_jansky(self) -> None:
        """
        Converts the y data and its errors from Watts per square meter to Jansky.
        """
        self._ydata = utils.watt_to_jansky(self._ydata, self._xdata)
        self._yerr = utils.watt_to_jansky(
            self._yerr, self._xdata) if self._yerr is not None else None

    def convert_to_watt(self) -> None:
        """
        Converts the y data and its errors from Jansky to Watts per square meter.
        """
        self._ydata = utils.jansky_to_watt(self._ydata, self._xdata)
        self._yerr = utils.jansky_to_watt(
            self._yerr, self._xdata) if self._yerr is not None else None

    def querry_vizier_data(self, radius: float = 5, target: str = None):
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
        table = utils.querry_vizier_data(radius, target)
        table.sort('sed_freq', reverse=True)
        return table
    
    def get_common_filters(self, table = None) -> dict:
        """
        Returns the common filters between the data and the bandpasses.

        Parameters:
        filter (array-like): The filters of the data.

        Returns:
        dict: The common filters.
        """
        bandpass_name = utils.get_bandpass_name()
        if table is None:
            table = self._table
            filter = table['sed_filter'].data
            return utils.get_common_filters(filter, bandpass_name)
        
        else : 
            table = unique(table, keys='sed_freq')
            table.sort('sed_freq', reverse=True)
            filter = table['sed_filter'].data
            return utils.get_common_filters(filter, bandpass_name)
    
    def restrict_data_vizier(self, table):
        """
        Restricts the data from a Vizier query based on the common bandpasses.

        Parameters:
        table (array-like): The data from the Vizier query.
        bandpass_name (list[str]): The names of the common bandpasses.
        """
        bandpass_name = utils.get_bandpass_name()
        filter = np.unique(table['sed_filter'].data)
        common_filters = utils.get_common_filters(filter, bandpass_name)
        mask = np.isin(table['sed_filter'].data, list(common_filters.keys()))
        table = table[mask]
        return table

    def set_vizier_data(self, table) -> None:
        """
        Sets the data from a Vizier query. Fill nan with 0 for compatibility with the rest of the code.
        The data is then aggregated by filter and the x and y data are set.

        Parameters:
        table (array-like): The data from the Vizier query.
        """
        data = np.array(table['sed_eflux'])
        has_nan = np.isnan(data)
        table[has_nan]['sed_eflux'] = 0.0

        table = utils.aggregate_table(table, column='sed_filter', fct=np.mean)

        self._table = table
        self._xdata = (np.asarray(table['sed_freq'])*table['sed_freq'].unit).to(u.um, equivalencies=u.spectral()).value
        self._ydata = table['sed_flux'].value
        self._yerr = table['sed_eflux'].value

    def restrict_data(self, ListOfCondition=list[str]):
        """
        Restricts the data based on a condition.
        Parameters:
        condition (str): The condition to restrict the data with.

        Examples:
        >>> data = Data(xdata=np.array([1, 2, 3, 4]), ydata=np.array([10, 20, 30, 40]))
        >>> data.restrict_data(['xdata > 2'])
        >>> print(data.get_xdata())
        [3 4]
        >>> print(data.get_ydata())
        [30 40]
        """
        for condition in ListOfCondition:
            condition = 'self._' + \
                condition.split(' ')[0] + ' ' + \
                ' '.join(condition.split(' ')[1:])
            restriction = eval(condition)

            self._xdata = self._xdata[restriction]
            self._ydata = self._ydata[restriction]
            self._xerr = self._xerr[restriction] if self._xerr is not None else None
            self._yerr = self._yerr[restriction] if self._yerr is not None else None
            
            if 'yerr' in condition:
                self._table = self._table[eval(condition.replace('yerr', "table['sed_eflux']"))] if self._table is not None else None
            if 'ydata' in condition:
                self._table = self._table[eval(condition.replace('ydata', "table['sed_flux']"))] if self._table is not None else None
            if 'xdata' in condition:
                condition = condition.replace('xdata', "table['sed_freq']")
                if '>' in condition and '=' not in condition:
                    condition = condition.replace('>', '<')
                elif '<' in condition and '=' not in condition:
                    condition = condition.replace('<', '>')
                elif '>=' in condition:
                    condition = condition.replace('>=', '<=')
                elif '<=' in condition:
                    condition = condition.replace('<=', '>=')
                limit = condition.split('>')[-1] if '>' in condition else condition.split('<')[-1]
                limit_value = (float(limit.strip('=<> ')) * u.um).to(u.GHz, equivalencies=u.spectral()).value
                condition = condition.replace(limit, str(limit_value))
                self._table = self._table[eval(condition)] if self._table is not None else None



    def add_data(self, xdata: np.array, ydata: np.array, xerr: np.array = None, yerr: np.array = None) -> None:
        """
        Adds new data to the existing data.

        Parameters:
        xdata (array-like): The new x data values.
        ydata (array-like): The new y data values.
        xerr (array-like, optional): The new errors in the x data values. Defaults to None.
        yerr (array-like, optional): The new errors in the y data values. Defaults to None.
        """
        if self._xdata is None:
            self._xdata = xdata
        else:
            self._xdata = np.concatenate((self._xdata, xdata))

        if self._ydata is None:
            self._ydata = ydata
        else:
            self._ydata = np.concatenate((self._ydata, ydata))

        if xerr is not None:
            if self._xerr is None:
                self._xerr = xerr
            else:
                self._xerr = np.concatenate((self._xerr, xerr))

        if yerr is not None:
            if self._yerr is None:
                self._yerr = yerr
            else:
                self._yerr = np.concatenate((self._yerr, yerr))

    def write_table_to_latex(self, Path: str, columns: list = None, column_names: list = None, wavelength: bool = True) -> None:
        """
        Writes the data table to a latex file.

        Parameters:
        Path (str): The file path to write the table to.
        """
        utils.write_table_to_latex(self._table, Path, columns, column_names, wavelength)

    def __str__(self):

        # if self._table is not None:
        #     old_stdout = sys.stdout
        #     sys.stdout = mystdout = StringIO()
        #     self._table[['_RAJ2000','_DEJ2000', '_tabname', 'sed_freq', 'sed_flux', 'sed_eflux', 'sed_filter']].pprint(max_lines=-1, max_width=-1)
        #     sys.stdout = old_stdout
        #     return mystdout.getvalue()
        table_str = "Data Table:\n"
        table_str += f"{'xdata':<20} {'ydata':<20} {'xerr':<20} {'yerr':<20} {'filter':<20}\n"
        table_str += "-"*90 + "\n"
        for i, x in enumerate(self._xdata):
            table_str += f"{format(x,'.3f'):<20} {format(self._ydata[i], '.3e'):<20} {format(self._xerr[i], '.3e') if self._xerr is not None else 'None':<20} {format(self._yerr[i], '.3e') if self._yerr is not None else 'None':<20} {str(self._table['sed_filter'][i]) if self._table is not None else 'None'}\n"
        return table_str
