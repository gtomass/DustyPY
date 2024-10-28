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

    def __init__(self, Model: Model = None , PATH: str = None, Lestimation: float = 1e4) -> None:
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
        self._DustyReconizer = {'Spectral': 'Spectral shape', # new black_body or engelke-marengo
                                'BB': 'Number of BB', # only with black_body
                                'Temperature': 'Temperature', 
                                'Luminosities': 'Luminosities', # only with black_body and BB > 1
                                'Absorption': 'SiO absorption depth', # new (not compatible with luminosities and only 1 star)
                                'Optical properties': 'optical properties index', # new
                                'Composition': 'Number of additional components', 
                                'Abundances': 'Abundances for these components',
                                'Size Distribution': 'SIZE DISTRIBUTION', # new MRN, Modified MRN
                                'Dust size': 'a(min) =',
                                'Sublimation temperature': 'Tsub', # new
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

    def AvailableComposition(self) -> list:
        """
        Returns a list of available compositions for the Dusty model.

        Returns:
        list: A list of available compositions.
        """
        return [
            file.split('/')[-1].split('.')[0]
            for file in glob.glob(f'{os.path.join(self._dustyPath,'data','Lib_nk','*.nk')}')
        ]
        

    def ChangeParameter(self) -> None:
        """
        Changes the parameters of the Dusty model based on the current model settings.
        """

        change = utils.build_change_dict(self._Model)
        utils.ChangeParameter(os.path.join(self._dustyPath,
                                           self._Model.get_Name(),
                                           self._Model.get_Name()+'.inp'),
                                           change=change,car=self._DustyReconizer,nstar=int(self._Model.get_NbStar()))

    def PrintParam(self) -> None:
        """
        Prints the current parameters of the Dusty model.
        """

        name = self._Model.get_Name()
        utils.PrintFile(os.path.join(self._dustyPath,name,name+'.inp'),stop=73)

    def LunchDusty(self,verbose: str = 2) -> None:
        """
        Runs the Dusty simulation with the current model settings.
        """
        if verbose not in [0,1,2]:
            raise ValueError('The verbose parameter must be 0, 1 or 2')
        subprocess.check_call([f'./dusty model.mas {verbose}'],cwd=self._dustyPath,shell=True)


    def GetResults(self) -> dict:
        """
        Retrieves the results of the Dusty simulation.

        Returns:
        dict: A dictionary containing the results of the Dusty simulation.
        """
         
        result_file = utils.LoadFile(os.path.join(self._dustyPath,self._Model.get_Name(),self._Model.get_Name()+'.out'))
        line = result_file[utils.SearchLine(result_file,'tau0')].split(' ')
        keys = utils.SuppCarList(line, ['###','','\n'])
        line = result_file[utils.SearchLine(result_file,'tau0') + 3].split(' ')
        values = [float(el) for el in utils.SuppCarList(line,['','1','\n'])]
        return utils.ListToDict(keys,values)
    
    def MakeSED(self,distance, Jansky: bool = True , um: bool = True) -> None:
        
        """
        Generates the Spectral Energy Distribution (SED) for the given distance.
        
        Parameters:
        distance (float): The distance in parsec to the object for which the SED is being generated.
        Jansky (bool, optional): If True, converts the flux to Jansky. Defaults to True.
        um (bool, optional): If True, keeps the wavelength in micrometers. If False, converts to meters. Defaults to True.
        """

        distance = distance * constants.pc

        results = self.GetResults()

        r_vrai = utils.CalculRayonVrai(results,self._Lest)
        FTot = utils.CalculFluxTotal(results['Fi(W/m2)'],r_vrai,distance)

        SED_file = utils.LoadFile(os.path.join(self._dustyPath,self._Model.get_Name(),self._Model.get_Name()+'.stb'))
        Wavelengths = utils.GetColumnSpectrum(SED_file,index = 0, index_header=6)
        Flux = FTot*utils.GetColumnSpectrum(SED_file,index=1,index_header=6)

        if Jansky:
            Flux = utils.WattToJansky(Flux,Wavelengths)
        if not um:
            Wavelengths = utils.UmToMeter(Wavelengths)

        self._SED.set_Wavelength(wavelength=Wavelengths)
        self._SED.set_Flux(Flux=Flux)

    def GetSED(self) -> SED:
        """
        Returns the Spectral Energy Distribution (SED) of the Dusty model.

        Returns:
        SED: The SED of the Dusty model.
        """

        return self._SED
    

    def MakeWavelength(self, number_of_wavelength: int = 200) -> None:
        """
        Generates a list of wavelengths for the Dusty model.

        Parameters:
        number_of_wavelength (int, optional): The number of wavelengths to generate. Defaults to 200.
        """

        wavelengths = utils.LogSpace(-2,4,number_of_wavelength)
        check = [wavelengths[i+1]/wavelengths[i] < 1.5 for i in range(len(wavelengths)-1)]
        if all(check):
            utils.WriteWavelength(os.path.join(self._dustyPath,'data','lambda_grid.dat'),wavelengths)
        else:
            raise ValueError('Not all wavelengths are spaced by less than 50%')

    
    def PlotSED(self, unit: str = None, xlim: tuple = None, ylim: tuple = None, ax = None, scale: str = 'linear', kwargs: dict = None) -> None:
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

        self._SED.PlotSED(unit=unit,xlim=xlim,ylim=ylim,ax=ax,scale=scale,kwargs=kwargs)
    

    def __Check(self) -> None:
        """
        Checks the consistency of the Dusty model's attributes.
        """
         
        for species in self._Model.get_Dust().get_Composition():
            if species not in self.AvailableComposition():
                raise ValueError(f'The following species does not exist: {species}')

    def __CreateDustyFile(self) -> None:
        """
        Creates the necessary Dusty model files based on the current model settings.
        """
        os.makedirs(os.path.join(self._dustyPath,self._Model.get_Name()), exist_ok=True)
        #print(os.path.join(os.path.dirname(__file__),'Mod.inp'))

        subprocess.call(
            [   
                f"cp {os.path.join(os.path.dirname(__file__),'Mod.inp')} {os.path.join(self._dustyPath ,self._Model.get_Name(), self._Model.get_Name()+'.inp')}"
            ],
            shell=True
        )

        with open(os.path.join(self._dustyPath,'model.mas'),'w') as file:
            file.write(f'{os.path.join(self._Model.get_Name(),self._Model.get_Name())+'.inp'}\n')