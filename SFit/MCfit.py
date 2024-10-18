import pymcmcstat
import pymcmcstat.MCMC
import utils
import Data
import numpy as np


class Fit():

    def __init__(self, Model = None, ParamFit = {}, Param={}):
        self._Model = Model
        self._ParamFit = ParamFit
        self._Param = Param

    def set_Model(self, Model):
        self._Model = Model

    def set_ParamFit(self, ParamFit):
        self._ParamFit = ParamFit

    def set_Param(self, Param):
        self._Param = Param

    def Fit(self, data = Data.Data):
        y= data.get_ydata()
        x= data.get_xdata()
        #sig = DF
        
        # x=x[x<10]
        # y=y[:len(x)]
        
        nds=len(x)
        
        results = {}
        mcstat = pymcmcstat.MCMC.MCMC()
        mcstat.data.add_data_set(x,y)
        utils.SetMCMCParam(mcstat,self._Param)
        mcstat.model_settings.define_model_settings(sos_function=utils.Chi2)
        mcstat.simulation_options.define_simulation_options(*self._ParamFit)

        mcstat.run_simulation()
        results = mcstat.simulation_results.results.copy()

        result = results
        chain = result['chain']
        s2chain = result['s2chain']
        sschain = result['sschain']
        names = result['names']
        # define burnin
        burnin = int(result['nsimu']/2)
        # display chain statistics
        mcstat.chainstats(chain[burnin:,:], result)
        return None

    
if __name__=='__main__':
    dat = Data.Data(xdata=[0,1,2,3,4,5,6,7,8,9,10],ydata=[0,1,2,3,4,5,6,7,8,9,10])

    ParamFit= {'nsimu':int(10000),
        'updatesigma':True,
        'method':'dram',
        'adaptint':100,
        'verbosity':0,
        'waitbar':True}

    f = Fit(Param={'I':{'theta0':1,'minimum':0,'maximum':2}})

    f.Fit(data=dat)
    