import SFit.stars as stars
import SFit.dusty as dusty
import SFit.utils as utils
import SFit.DustyFit as DustyFit
import SFit.Data as Data
import matplotlib.pyplot as plt

import numpy as np

if __name__== '__main__':

    distance = 1000*3.26*365*86400*3e8
    c=3e8

    #print(7.6/(7.6+9.8))

    S1 = stars.Star(Name='E1',Temperature=2441.62,Luminosity=0.37)
    S2 = stars.Star(Name='E2',Temperature=8072.59,Luminosity=0.63)

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

    dustyMod.ChangeParameter()
    dustyMod.LunchDusty()

    kwargs_data = {'marker':'+','color':'r'}
    kwargs_dusty = {'color':'k'}

    Dat = Data.Data()
    d = Dat.ImportData('/Users/gabriel/Documents/Stage/data/spectre.csv').fillna(0)

    d = d[1e-26*d._sed_freq*1e9*d._sed_flux>1e-12]
    d = d[1e-26*d._sed_freq*1e9*d._sed_flux<1e-10]
    d = d[c/(d._sed_freq*1e9) < 10e-6]
    
    fig,ax = plt.subplots()

    Dat.set_xdata(np.asarray(3e8/(d['_sed_freq']*1e9)*1e6))
    Dat.set_ydata(np.asarray(d['_sed_flux']))
    #print(Dat.get_ydata())
    Dat.ConvertWatt()
    Dat.UnredData(EBV=1.07)
    Dat.ConvertJansky()
    # Dat.ScatterData(xlim=(0,5),ylim=(0,100), ax=ax, kwargs=kwargs_data)
    # dustyMod.MakeSED(distance=distance)
    # dustyMod.PlotSED(xlim=(0,5),ylim=(0,100), ax=ax, kwargs=kwargs_dusty)
    # plt.show()


    fit = DustyFit.DustyFit(dustyMod, Data=Dat)

    Param = {'Temp1':{'theta0':3000,'minimum':2000,'maximum':4000},
             'Lum1':{'theta0': 0.5, 'minimum':0,'maximum':1},
             'Temp2':{'theta0':7350,'minimum':5000,'maximum':10000},
             'Lum2':{'theta0':0.5,'minimum':0,'maximum':1},
             'Ltot':{'theta0':1,'minimum':0,'maximum':10},
            }
    
    ParamFit = {
                        'nsimu': 100,
                        'updatesigma': True,
                        'method': 'dram',
                        'adaptint': 100,
                        'verbosity': 0,
                        'waitbar': True,
                    }
    fit.set_Param(Param)
    fit.set_ParamFit(ParamFit)

    fit.LunchFit()
    fit.get_Fit().PrintResults()
    #fit.PlotResults(xlim=(0,5),ylim=(0,100), ax=ax, kwargs_fit=kwargs_dusty, kwargs_data=kwargs_data)



