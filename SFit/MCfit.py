import pymcmcstat
import pymcmcstat.MCMC
import SFit.utils as utils
import SFit.Data as Data
import numpy as np


class Fit():

    def __init__(self, Data = Data.Data, Model = pymcmcstat.MCMC.MCMC(), ParamFit = {}, Param={}):
        self._Model = Model
        self._Data = Data
        self._ParamFit = ParamFit
        self._Param = Param
        self._Results = {}

    def set_Data(self, Data = Data.Data):
        self._Data = Data

    def get_Data(self):
        return self._Data

    def set_Model(self):
        utils.SetMCMCParam(self._Model,self._Param)
        self._Model.simulation_options.define_simulation_options(*self._ParamFit)

    def set_ParamFit(self, ParamFit):
        self._ParamFit = ParamFit
        
    def set_Param(self, Param):
        self._Param = Param

    def set_Chi2Func(self, func):
        self._Model.model_settings.define_model_settings(sos_function=func)

    def get_Results(self):
        return self._Results

    def Fit(self, Chi2):
        y= self._Data.get_ydata()
        x= self._Data.get_xdata()
        self._Model.data.add_data_set(x,y)
        self.set_Model()

        self.set_Chi2Func(Chi2)
        self._Model.run_simulation()

        results = {}
        results = self._Model.simulation_results.results.copy()

        self._Results = results

    def PrintResults(self):
        chain = self._Results['chain']
        burnin = int(self._Results['nsimu']/2)
        self._Model.chainstats(chain[burnin:,:], self._Results)

    def PredictionModel(self):
        nds = len(self._Data.get_xdata())
        print(nds)
        def pred_modelfun(preddata, theta):
            return utils.model(theta[:2], preddata.xdata[0],preddata.ydata[0]).reshape(nds,)
        self._Model.PI.setup_prediction_interval_calculation(
        results=self._Results,
        data=self._Model.data,
        modelfunction=pred_modelfun)
        self._Model.PI.generate_prediction_intervals()
        self._Model.PI.plot_prediction_intervals(adddata=True, figsizeinches=(6, 6))

    def get_ModelSpline():
        pass 
    

if __name__=='__main__':
    dat = Data.Data(xdata=[0,1,2,3,4,5,6,7,8,9,10],ydata=[0,1,2,3,4,5,6,7,8,9,10])

    ParamFit= {'nsimu':int(10000),
        'updatesigma':True,
        'method':'dram',
        'adaptint':100,
        'verbosity':0,
        'waitbar':True}

    f = Fit(Param={'I':{'theta0':1,'minimum':0,'maximum':2}},Data=dat)

    f.Fit()
    f.PrintResults()
    result = f.get_Results()['theta'][0]
    x = np.linspace(0,10,100)
    y = utils.model(q=result,x=dat.get_xdata(),y=dat.get_ydata())

    utils.Plot(y,dat.get_xdata())

    