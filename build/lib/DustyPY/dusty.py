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

    def __init__(self, PATH='', Model=Model(), Lestimation=1e4):
        """
        Initializes an instance of the Dusty class.

        Parameters:
        PATH (str, optional): The path to the Dusty model files. Defaults to an empty string.
        Model (Model, optional): The model used in the Dusty simulation. Defaults to an instance of Model.
        Lestimation (float, optional): The luminosity estimation. Defaults to 1e4.
        """
        self._dustyPath = PATH
        self._Model = Model
        self._Lest = Lestimation
        self._DustyReconizer = {'BB': 'Number of BB',
                                'Temperature': 'Temperature', 
                                'Luminosities': 'Luminosities',
                                'Opacity': 'tau(min)',
                                'Composition': 'Number of additional components',
                                'Abundances': 'Abundances for these components'}
        self._SED = SED.SED()
        self.__Check()
        self.__CreateDustyFile()

    def set_Model(self, Model):
        """
        Sets the model used in the Dusty simulation.

        Parameters:
        Model (Model): The model to be used in the Dusty simulation.
        """
        self._Model = Model
        self.__Check()

    def get_Model(self):
        """
        Returns the model used in the Dusty simulation.

        Returns:
        Model: The model used in the Dusty simulation.
        """
        return self._Model

    def set_PATH(self, Path):
        """
        Sets the path to the Dusty model files.

        Parameters:
        Path (str): The path to the Dusty model files.
        """
        self._dustyPath = Path

    def AvailableComposition(self):
        """
        Returns a list of available compositions for the Dusty model.

        Returns:
        list: A list of available compositions.
        """
        return [
            file.split('/')[-1].split('.')[0]
            for file in glob.glob(f'{self._dustyPath}/Lib_nk/*.nk')
        ]
        

    def ChangeParameter(self):
        """
        Changes the parameters of the Dusty model based on the current model settings.
        """

        name = self._Model.get_Name()
        Stars = self._Model.get_Stars()

        T = [str(Star.get_Temperature()) for Star in Stars]
        L = [str(Star.get_Luminosity())for Star in Stars]

        dust = self._Model.get_Dust()
        comp = "\n        ".join(
            f"{f'Lib_nk/{comp}'}.nk" for comp in dust.get_Composition().keys()
        )
        nbcomp = str(len(dust.get_Composition().keys()))
        abondances = ", ".join(f'{ab}' for ab in dust.get_Composition().values())

        change = {  'BB': f'        	Number of BB = {len(T)} \n',
                    'Temperature': f'        	Temperature = {', '.join(T)} K \n', 
                    'Luminosities': f'        	Luminosities = {', '.join(L)} \n',
                    'Opacity': f'        - tau(min) = {dust.get_tau()}; tau(max) = {dust.get_tau()}  % for the visual wavelength \n' ,
                    'Composition': f'	Number of additional components = {nbcomp} properties listed in: \n        {comp}\n',
                    'Abundances': f'   Abundances for these components = {abondances} \n'    ,
                }

        utils.ChangeParameter(self._dustyPath+name+'.inp',change=change,car=self._DustyReconizer,ncomp=int(nbcomp))

    def PrintParam(self):
        """
        Prints the current parameters of the Dusty model.
        """

        name = self._Model.get_Name()
        utils.PrintFile(self._dustyPath+name+'.inp',stop=73)

    def LunchDusty(self):
        """
        Runs the Dusty simulation with the current model settings.
        """

        subprocess.check_call(['./dusty'],cwd=self._dustyPath)


    def GetResults(self):
        """
        Retrieves the results of the Dusty simulation.

        Returns:
        dict: A dictionary containing the results of the Dusty simulation.
        """
         
        result_file = utils.LoadFile(self._dustyPath+self._Model.get_Name()+'.out')
        line = result_file[utils.SearchLine(result_file,'tau0')].split(' ')
        keys = utils.SuppCarList(line, ['###','','\n'])
        line = result_file[utils.SearchLine(result_file,'tau0') + 3].split(' ')
        values = [float(el) for el in utils.SuppCarList(line,['','1','\n'])]
        return utils.ListToDict(keys,values)
    
    def MakeSED(self,distance, Jansky = True , um = True):
        
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
        FTot = utils.CalculFluxTotal(results['F1(W/m2)'],r_vrai,distance)

        SED_file = utils.LoadFile(self._dustyPath+self._Model.get_Name()+'.stb')
        Wavelengths = utils.GetColumnSpectrum(SED_file,index = 0, index_header=4)
        Flux = FTot*utils.GetColumnSpectrum(SED_file,index=1,index_header=4)

        if Jansky:
            Flux = utils.WattToJansky(Flux,Wavelengths)
        if not um:
            Wavelengths = utils.UmToMeter(Wavelengths)

        self._SED.set_Wavelength(wavelength=Wavelengths)
        self._SED.set_Flux(Flux=Flux)

    def GetSED(self):
        """
        Returns the Spectral Energy Distribution (SED) of the Dusty model.

        Returns:
        SED: The SED of the Dusty model.
        """

        return self._SED
    
    def PlotSED(self, unit=None, xlim=None, ylim=None, ax=None, scale='linear', kwargs=None):
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
    

    def __Check(self):
        """
        Checks the consistency of the Dusty model's attributes.
        """
         
        for species in self._Model.get_Dust().get_Composition():
            if species not in self.AvailableComposition():
                raise ValueError(f'The following species does not exist: {species}')

    def __CreateDustyFile(self):
        """
        Creates the necessary Dusty model files based on the current model settings.
        """

        subprocess.call(
            [
                f"cp {os.getcwd()}/SFit/Mod.inp {self._dustyPath + self._Model.get_Name() + '.inp'}"
            ],
            shell=True,
        )