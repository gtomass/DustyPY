import DustyPY.stars as stars
import DustyPY.dusty as dusty
import DustyPY.utils as utils
import DustyPY.DustyFit as DustyFit
import DustyPY.Data as Data
import matplotlib.pyplot as plt

import numpy as np

if __name__== '__main__':

    distance = 1000

    #print(7.6/(7.6+9.8))

    S1 = stars.Star(Name='E1',Temperature=3000,Luminosity=9.8/(7.6+9.8))
    S2 = stars.Star(Name='E2',Temperature=7250,Luminosity=7.6/(7.6+9.8))

    Composition = {'Mg2SiO4':0.04, 
                   'MgSiO3':0.1, 
                   'H2Oice': 0.05,
                   'Al2O3-comp': 0.15,
                   'MgFeSiO4':0.65,
                   'FeO':0.01 
                   }

    dust = stars.Dust(DustSize={'amin' : 0.005, 'amax': 0.25}, Composition=Composition,tau=0)

    mod = stars.Model('AFGL4106', NbStar=2, Stars=[S1,S2], Dust=dust, distance=distance)
   
    dustyMod = dusty.Dusty(PATH='/Users/gabriel/Documents/Stage/code/dustyV2/', Model=mod)

    #dustyMod.ChangeParameter()
    #dustyMod.LunchDusty()

    kwargs_data = {'fmt':'+','color':'r'}
    kwargs_dusty = {'color':'k'}

    Dat = Data.Data()

    table = Dat.QuerryVizierData(radius = 5, target='HD302821')
    Dat.SetVizierData(table)
    Dat.RestrictData(['yerr != 0.', 'xdata < 10'])
    
    fig,ax = plt.subplots()
    Dat.ConvertWatt()
    Dat.UnredData(EBV=1.07)
    Dat.ConvertJansky()
    Dat.RestrictData(['ydata < 100'])
    Dat.ScatterData(xlim=(0,5),ylim=(0,100), ax=ax, kwargs=kwargs_data)
    dustyMod.MakeSED(distance=distance)
    dustyMod.PlotSED(xlim=(0,5),ylim=(0,100), ax=ax, kwargs=kwargs_dusty)
    plt.show()


    # fit = DustyFit.DustyFit(dustyMod, Data=Dat)

    # Param = {'Temp1':{'theta0':3000,'minimum':2000,'maximum':4000},
    #          'Lum1':{'theta0': 0.5, 'minimum':0,'maximum':1},
    #          'Temp2':{'theta0':7350,'minimum':5000,'maximum':10000},
    #          'Lum2':{'theta0':0.5,'minimum':0,'maximum':1},
    #          'Ltot':{'theta0':1,'minimum':0,'maximum':10},
    #         }
    
    # ParamFit = {
    #                     'nsimu': 100,
    #                     'updatesigma': True,
    #                     'method': 'dram',
    #                     'adaptint': 100,
    #                     'verbosity': 0,
    #                     'waitbar': True,
    #                 }
    # fit.set_Param(Param)
    # fit.set_ParamFit(ParamFit)

    # fit.LunchFit()
    # fit.get_Fit().PrintResults()
    # fit.PlotResults(xlim=(0,5),ylim=(0,100), ax=ax, kwargs_fit=kwargs_dusty, kwargs_data=kwargs_data)



