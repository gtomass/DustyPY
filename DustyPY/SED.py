from . import utils as utils


class SED():
    """
    Class representing a Spectral Energy Distribution (SED).

    Attributes:
    _Flux (list): The fluxes of the SED.
    _Wavelength (list): The wavelengths of the SED.
    """

    def __init__(self, Flux=list, wavelength=list) -> None:
        """
        Initializes an instance of the SED class.

        Parameters:
        Flux (list, optional): The fluxes of the SED. Defaults to an empty list.
        wavelength (list, optional): The wavelengths of the SED. Defaults to an empty list.
        """
        self._Flux = Flux
        self._Wavelength = wavelength

    def set_Flux(self, Flux: list) -> None:
        """
        Sets the fluxes of the SED.

        Parameters:
        Flux (list): The fluxes of the SED.
        """
        self._Flux = Flux

    def set_Wavelength(self, wavelength: list) -> None:
        """
        Sets the wavelengths of the SED.

        Parameters:
        wavelength (list): The wavelengths of the SED.
        """
        self._Wavelength = wavelength

    def get_Flux(self) -> list:
        """
        Returns the fluxes of the SED.

        Returns:
        list: The fluxes of the SED.
        """
        return self._Flux

    def get_Wavelength(self) -> list:
        """
        Returns the wavelengths of the SED.

        Returns:
        list: The wavelengths of the SED.
        """
        return self._Wavelength
    
    def integrate_bandpass(self, bandpass: list) -> float:
        """
        Integrates the SED over a bandpass.

        Parameters:
        bandpass (list): The
        """
        return utils.integrate_SED_bandpass(self._Wavelength, self._Flux, bandpass)

    def plot_SED(self, unit: str = None, xlim: tuple = None, ylim: tuple = None, ax=None, scale: str = 'linear', kwargs: dict = None, normalize: bool = False) -> None:
        """
        plots the SED using the fluxes and wavelengths.

        Parameters:
        unit (str, optional): The unit of the axes. Defaults to None.
        xlim (tuple, optional): The limits of the x-axis. Defaults to None.
        ylim (tuple, optional): The limits of the y-axis. Defaults to None.
        ax (matplotlib.axes.Axes, optional): The axis on which to plot. Defaults to None.
        scale (str, optional): The scale of the axes ('linear' or 'log'). Defaults to 'linear'.
        kwargs (dict, optional): Additional arguments for the plot function. Defaults to None.
        """
        utils.plot(self._Flux, self._Wavelength, unit=unit, xlim=xlim,
                   ylim=ylim, ax=ax, scale=scale, kwargs=kwargs, normalize=normalize)
        
    def scatter_SED_bandpass(self,bandpass,unit: str = None, xlim: tuple = None, ylim: tuple = None, ax=None, scale: str = 'linear', kwargs: dict = None, normalize: bool = False) -> None:
        """
        Scatter the SED over a bandpass.

        Parameters:
        bandpass (list): The
        ax (matplotlib.axes.Axes, optional): The axis on which to plot. Defaults to None.
        kwargs (dict, optional): Additional arguments for the scatter function. Defaults to None.
        """
        central =[utils.get_central_wavelegnth(utils.get_bandpass(f))/10000 for f in bandpass.values()]
        f = self.integrate_bandpass(bandpass)
        utils.scatter_plot(self.integrate_bandpass(bandpass),central, unit=unit, xlim=xlim,
                   ylim=ylim, ax=ax, scale=scale, kwargs=kwargs, normalize=normalize)

