import subprocess
import glob
import SFit.utils as utils 
from SFit.stars import Model
import os



class Dusty():

    def __init__(self, PATH='', Model=Model):
        self._dustyPath = PATH
        self._Model = Model
        self._DustyReconizer = {'BB': 'Number of BB',
                                'Temperature': 'Temperature', 
                                'Luminosities': 'Luminosities',
                                'Opacity': 'tau(min)',
                                'Composition': 'Number of additional components',
                                'Abundances': 'Abundances for these components'}
        self.__Check()
        self.__CreateDustyFile()

    def set_Model(self, Model):
        self._Model = Model
        self.__Check()

    def get_Model(self):
        return self._Model

    def set_PATH(self, Path):
        self._dustyPath = Path

    def AvailableComposition(self):
        NkFiles = [file.split('/')[-1].split('.')[0] for file in glob.glob(self._dustyPath+'/Lib_nk/*.nk')]
        return NkFiles
        

    def ChangeParameter(self):

        name = self._Model.get_Name()
        Stars = self._Model.get_Stars()

        T = [str(Star.get_Temperature()) for Star in Stars]
        L = [str(Star.get_Luminosity())for Star in Stars]

        dust = self._Model.get_Dust()
        comp = " \n        ".join(f'{'Lib_nk/'+comp}'+'.nk' for comp in dust.get_Composition().keys())
        nbcomp = str(len(dust.get_Composition().keys()))
        abondances = ", ".join(f'{ab}' for ab in dust.get_Composition().values())
        
        change = {  'BB': f'        	Number of BB = {len(T)} \n',
                    'Temperature': f'        	Temperature = {', '.join(T)} K \n', 
                    'Luminosities': f'        	Luminosities = {', '.join(L)} \n',
                    'Opacity': f'        - tau(min) = {dust.get_tau()}; tau(max) = {dust.get_tau()}  % for the visual wavelength \n' ,
                    'Composition': f'	Number of additional components = {nbcomp} properties listed in: \n        {comp}',
                    'Abundances': f'   Abundances for these components = {abondances} \n \n'    ,
                }

        utils.ChangeParameter(self._dustyPath+name+'.inp',change=change,car=self._DustyReconizer)

    def PrintParam(self):
        name = self._Model.get_Name()
        utils.PrintFile(self._dustyPath+name+'.inp',stop=73)

    def LunchDusty(self):
        subprocess.check_call(['./dusty'],cwd=self._dustyPath)

    def GetResults(self):
        print('Not implemented yet')
        return None

    def __Check(self):
        for species in self._Model.get_Dust().get_Composition():
            if species not in self.AvailableComposition():
                raise ValueError(f'The following species does not exist: {species}')

        if sum([Star.get_Luminosity() for Star in self._Model.get_Stars()]) != 1.:
            raise Exception('Sum of Luminosities must be 1')

    def __CreateDustyFile(self):
        if self._dustyPath+self._Model.get_Name()+'.inp' in glob.glob(self._dustyPath):
            pass
        else:
            subprocess.call([f'cp {os.getcwd()+'/SFit/Mod.inp'} {self._dustyPath+self._Model.get_Name()+'.inp'}'],shell=True)