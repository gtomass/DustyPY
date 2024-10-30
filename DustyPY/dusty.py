import subprocess
import glob
from . import utils as utils
from .stars import Model as Model
from . import SED as SED
from . import constants
import os


class Dusty():
    """
    Class representing a Dusty model with a specified path, model, and luminosity estimation.

    Attributes:
    _dustyPath (str): The path to the Dusty model files.
    _Model (Model): The model used in the Dusty simulation.
    _Lest (float): The luminosity estimation of the model.
    _DustyReconizer (dict): A dictionary mapping Dusty parameters to their descriptions.
    _SED (SED): The Spectral Energy Distribution (SED) associated with the Dusty model.
    """

    def __init__(self, Model: Model = None, PATH: str = None, Lestimation: float = 1e4) -> None:
        """
        Initializes an instance of the Dusty class.

        Parameters:
        PATH (str, optional): The path to the Dusty model files. Defaults to an empty string.
        Model (Model, optional): The model used in the Dusty simulation. Defaults to an instance of Model.
        Lestimation (float, optional): The luminosity estimation. Defaults to 1e4.
        """

        if PATH is None:
            PATH = ''
        if Model is None:
            Model = Model()

        self._dustyPath = PATH
        self._Model = Model
        self._Lest = Lestimation
        self._DustyReconizer = {'Spectral': 'Spectral shape',
                                'BB': 'Number of BB',
                                'Temperature': 'Temperature',
                                'Luminosities': 'Luminosities',
                                'Absorption': 'SiO absorption depth',
                                'Optical properties': 'optical properties index',
                                'Composition': 'Number of additional components',
                                'Abundances': 'Abundances for these components',
                                'Size Distribution': 'SIZE DISTRIBUTION',
                                'Dust size': 'a(min) =',
                                'Sublimation temperature': 'Tsub',
                                'Opacity': 'tau(min)'
                                }
        self._SED = SED.SED()
        self.__Check()
        self.__CreateDustyFile()

    def set_Model(self, Model: Model) -> None:
        """
        Sets the model used in the Dusty simulation.

        Parameters:
        Model (Model): The model to be used in the Dusty simulation.
        """
        self._Model = Model
        self.__Check()

    def get_Model(self) -> Model:
        """
        Returns the model used in the Dusty simulation.

        Returns:
        Model: The model used in the Dusty simulation.
        """
        return self._Model

    def set_PATH(self, Path: str) -> None:
        """
        Sets the path to the Dusty model files.

        Parameters:
        Path (str): The path to the Dusty model files.
        """
        self._dustyPath = Path

    def get_available_composition(self) -> list:
        """
        Returns a list of available compositions for the Dusty model.

        Returns:
        list: A list of available compositions.
        """
        return [
            file.split('/')[-1].split('.')[0]
            for file in glob.glob(f'{os.path.join(self._dustyPath, 'data', 'Lib_nk', '*.nk')}')
        ]

    def change_parameter(self) -> None:
        """
        Changes the parameters of the Dusty model based on the current model settings.
        """

        change = utils.build_change_dict(self._Model)
        utils.change_parameter(os.path.join(self._dustyPath,
                                            self._Model.get_Name(),
                                            self._Model.get_Name()+'.inp'),
                               change=change, car=self._DustyReconizer, nstar=int(self._Model.get_NbStar()))

    def print_param(self) -> None:
        """
        Prints the current parameters of the Dusty model.
        """

        name = self._Model.get_Name()
        utils.print_file(os.path.join(
            self._dustyPath, name, name+'.inp'), stop=73)

    def lunch_dusty(self, verbose: str = 2) -> None:
        """
        Runs the Dusty simulation with the current model settings.
        """
        if verbose not in [0, 1, 2]:
            raise ValueError('The verbose parameter must be 0, 1 or 2')
        subprocess.check_call(
            [f'./dusty model.mas {verbose if verbose != None else ''}'], cwd=self._dustyPath, shell=True)

    def get_results(self) -> dict:
        """
        Retrieves the results of the Dusty simulation.

        Returns:
        dict: A dictionary containing the results of the Dusty simulation.
        """

        result_file = utils.load_file(os.path.join(
            self._dustyPath, self._Model.get_Name(), self._Model.get_Name()+'.out'))
        line = result_file[utils.search_line(result_file, 'tau0')].split(' ')
        keys = utils.supp_car_list(line, ['###', '', '\n'])
        line = result_file[utils.search_line(
            result_file, 'tau0') + 3].split(' ')
        values = [float(el)
                  for el in utils.supp_car_list(line, ['', '1', '\n'])]
        return utils.list_to_dict(keys, values)

    def make_SED(self, distance, Jansky: bool = True, um: bool = True) -> None:
        """
        Generates the Spectral Energy Distribution (SED) for the given distance.

        Parameters:
        distance (float): The distance in parsec to the object for which the SED is being generated.
        Jansky (bool, optional): If True, converts the flux to Jansky. Defaults to True.
        um (bool, optional): If True, keeps the wavelength in micrometers. If False, converts to meters. Defaults to True.
        """

        distance = distance * constants.pc

        results = self.get_results()

        r_vrai = utils.calcul_rayon_vrai(results, self._Lest)
        FTot = utils.calcul_flux_total(results['Fi(W/m2)'], r_vrai, distance)

        SED_file = utils.load_file(os.path.join(
            self._dustyPath, self._Model.get_Name(), self._Model.get_Name()+'.stb'))
        Wavelengths = utils.get_column_spectum(
            SED_file, index=0, index_header=6)
        Flux = FTot*utils.get_column_spectum(SED_file, index=1, index_header=6)

        if Jansky:
            Flux = utils.watt_to_jansky(Flux, Wavelengths)
        if not um:
            Wavelengths = utils.UmToMeter(Wavelengths)

        self._SED.set_Wavelength(wavelength=Wavelengths)
        self._SED.set_Flux(Flux=Flux)

    def get_SED(self) -> SED:
        """
        Returns the Spectral Energy Distribution (SED) of the Dusty model.

        Returns:
        SED: The SED of the Dusty model.
        """

        return self._SED

    def make_wavelength(self, number_of_wavelength: int = 200) -> None:
        """
        Generates a list of wavelengths for the Dusty model.

        Parameters:
        number_of_wavelength (int, optional): The number of wavelengths to generate. Defaults to 200.
        """

        wavelengths = utils.log_space(-2, 4, number_of_wavelength)
        check = [wavelengths[i+1]/wavelengths[i] <
                 1.5 for i in range(len(wavelengths)-1)]
        if all(check):
            utils.write_wavelength(os.path.join(
                self._dustyPath, 'data', 'lambda_grid.dat'), wavelengths)
        else:
            raise ValueError('Not all wavelengths are spaced by less than 50%')

    def plot_SED(self, unit: str = None, xlim: tuple = None, ylim: tuple = None, ax=None, scale: str = 'linear', kwargs: dict = None) -> None:
        """
        Plots the Spectral Energy Distribution (SED) of the Dusty model.

        Parameters:
        unit (str, optional): The unit of the axes. Defaults to None.
        xlim (tuple, optional): The limits of the x-axis. Defaults to None.
        ylim (tuple, optional): The limits of the y-axis. Defaults to None.
        ax (matplotlib.axes.Axes, optional): The axis on which to plot. Defaults to None.
        scale (str, optional): The scale of the axes ('linear' or 'log'). Defaults to 'linear'.
        kwargs (dict, optional): Additional arguments for the plot function. Defaults to None.
        """

        self._SED.plot_SED(unit=unit, xlim=xlim, ylim=ylim,
                           ax=ax, scale=scale, kwargs=kwargs)

    def __Check(self) -> None:
        """
        Checks the consistency of the Dusty model's attributes.
        """

        for species in self._Model.get_Dust().get_Composition():
            if species not in self.get_available_composition():
                raise ValueError(
                    f'The following species does not exist: {species}')

    def __CreateDustyFile(self) -> None:
        """
        Creates the necessary Dusty model files based on the current model settings.
        """
        os.makedirs(os.path.join(self._dustyPath,
                    self._Model.get_Name()), exist_ok=True)
        # print(os.path.join(os.path.dirname(__file__),'Mod.inp'))

        subprocess.call(
            [
                f"cp {os.path.join(os.path.dirname(__file__), 'Mod.inp')} {os.path.join(
                    self._dustyPath, self._Model.get_Name(), self._Model.get_Name()+'.inp')}"
            ],
            shell=True
        )

        with open(os.path.join(self._dustyPath, 'model.mas'), 'w') as file:
            file.write(f'{os.path.join(self._Model.get_Name(),
                       self._Model.get_Name())+'.inp'}\n')
