from . import dusty 
from . import MCfit
from . import Data
from . import utils as utils
import numpy as np
import matplotlib.pyplot as plt

class DustyFit():

    def __init__(self, Dusty, Data, ParamFit = None, Param = None, Fit = MCfit.Fit()):
        
        self._Dusty = Dusty
        self._Data = Data
        self._ParamFit = ParamFit
        self._Param = Param
        self._Fit = Fit

    def __InitFit(self):
        self._Fit.set_Data(Data=self._Data)
        self._Fit.set_ParamFit(ParamFit=self._ParamFit)
        self._Fit.set_Param(Param=self._Param)

    # Méthodes set
    def set_Dusty(self, Dusty):
        self._Dusty = Dusty

    def set_Data(self, Data):
        self._Data = Data

    def set_ParamFit(self, ParamFit):
        self._ParamFit = ParamFit

    def set_Param(self, Param):
        self._Param = Param

    def set_Fit(self, Fit):
        self._Fit = Fit

    # Méthodes get
    def get_Dusty(self):
        return self._Dusty

    def get_Data(self):
        return self._Data

    def get_ParamFit(self):
        return self._ParamFit

    def get_Param(self):
        return self._Param

    def get_Fit(self):
        return self._Fit
    
    def __setChange(self,change):
        for key in change.keys():
            if 'Temp' in key:
                self._Dusty.get_Model().get_Stars()[int(key.split('Temp')[-1])-1].set_Temperature(change[key])
            elif 'Lum' in key:
                self._Dusty.get_Model().get_Stars()[int(key.split('Lum')[-1])-1].set_Luminosity(change[key])
            elif key in ['Opacity']:
                self._Dusty.get_Model().get_Dust().set_Opacity(change[key])
            elif key in ['Composition','Abundances']:
               raise NotImplementedError('Composition and Abundances are not fittable yet')
    
    def __Chi2Dusty(self,theta,data):
        
        change = utils.ListToDict(list(self._Fit.get_Param().keys()),theta)

        self.__setChange(change)
        
        self._Dusty.ChangeParameter()
        self._Dusty.LunchDusty()
        self._Dusty.MakeSED(distance = self._Dusty.get_Model().get_Distance())

        self._Data.set_xdata(self._Dusty.GetSED().get_Wavelength())
        self._Data.set_ydata(self._Dusty.GetSED().get_Flux())


        try:
            xdata = data.xdata[0]
            ydata = data.ydata[0]
        except AttributeError:
            xdata,ydata = data.get_xdata(),data.get_ydata()

        ymodel = utils.model(theta[-1], xdata, self._Data.get_xdata(), self._Data.get_ydata()).reshape(ydata.shape)

        if self._Data.get_yerr() is not None:
            return np.nansum(((ymodel - ydata)/self._Data.get_yerr())**2)
        else:    
            return np.nansum((ymodel - ydata)**2)
            
    
    def LunchFit(self):
        self.__InitFit()
        self._Fit.Fit(Chi2=self.__Chi2Dusty)

    def PrintResults(self):
        self._Fit.PrintResults()
    
    def PlotResults(self,unit=None, xlim=None, ylim=None, ax=None, scale='linear', kwargs_fit=None, kwargs_data=None):
        result = self._Fit.get_Results()['theta']
        self.__setChange(utils.ListToDict(list(self._Fit.get_Param().keys()),result))
        print(utils.ListToDict(list(self._Fit.get_Param().keys()),result))
        self._Dusty.ChangeParameter()
        self._Dusty.LunchDusty()
        self._Dusty.MakeSED(distance = self._Dusty.get_Model().get_Distance())
        SED = self._Dusty.GetSED()
        utils.Plot(SED.get_Flux(),SED.get_Wavelength(),unit=unit,xlim=xlim,ylim=ylim,ax=ax,scale=scale,kwargs=kwargs_fit)
        if self._Data.get_yerr() is not None:
            utils.ErrorPlot(self._Data.get_ydata(),self._Data.get_xdata(),self._Data.get_yerr(),unit=unit,xlim=xlim,ylim=ylim,ax=ax,scale=scale,kwargs=kwargs_data)
        else:
            utils.ScatterPlot(self._Data.get_ydata(),self._Data.get_xdata(),unit=unit,xlim=xlim,ylim=ylim,ax=ax,scale=scale,kwargs=kwargs_data)
        plt.show()


    def __Check(self):
        pass
