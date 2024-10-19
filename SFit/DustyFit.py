import SFit.dusty
import SFit.MCfit
import SFit.Data
import SFit.utils as utils

class DustyFit():

    def __init__(self, Dusty, Data, ParamFit, Param, Fit = SFit.MCfit.Fit):
        
        self._Dusty = Dusty
        self._Data = Data
        self._ParamFit = ParamFit
        self._Param = Param
        self._Fit = Fit

    def InitFit(self):
        self._Fit.set_Data(Data=self._Data)
        self._Fit.set_ParamFit(ParamFit=self._ParamFit) 
        
    def InitParam(self):
        star = [star for star in self._Dusty.get_Model().get_Stars()]
        Param = utils.ListToDict()
        print(Param)
        self._Fit.set_Param(Param=self._Param)

    
    
