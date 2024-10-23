import DustyPY.stars as stars
import DustyPY.dusty as dusty
import DustyPY.utils as utils
import DustyPY.DustyFit as DustyFit
import DustyPY.Data as Data
import matplotlib.pyplot as plt

import numpy as np

if __name__== '__main__':

    distance = 197.0

    #print(7.6/(7.6+9.8))

    S1 = stars.Star(Name='E1',Temperature=3500,Luminosity=1)
    #S2 = stars.Star(Name='E2',Temperature=7250,Luminosity=7.6/(7.6+9.8))

    Composition = {'Al2O3-comp':0.2,'MgFeSiO4':0.16,'Ca2Al2SiO7':0.64}

    dust = stars.Dust(DustSize={'amin' : 0.01, 'amax': 1}, Composition=Composition,tau=0)

    mod = stars.Model('Prout', NbStar=1, Stars=[S1], Dust=dust, distance=distance)
   
    dustyMod = dusty.Dusty(PATH='/Users/gabriel/Documents/Stage/code/dustyV2/', Model=mod,Lestimation=1.8e5)
    print(dustyMod.AvailableComposition())

    dustyMod.ChangeParameter()
    dustyMod.LunchDusty()

    kwargs_data = {'fmt':'+','color':'r'}
    kwargs_dusty = {'color':'k'}

    Dat = Data.Data()

    table = Dat.QuerryVizierData(radius = 5, target='Betelgeuse')
    Dat.SetVizierData(table)
    Dat.RestrictData(['yerr != 0.'])
    
    fig,ax = plt.subplots()
    Dat.ConvertWatt()
    Dat.UnredData(EBV=0.4)
    Dat.ConvertJansky()
    # Dat.RestrictData(['ydata < 100'])
    Dat.ScatterData(xlim=(0,20),ylim=(0,6e4), ax=ax, kwargs=kwargs_data)
    dustyMod.MakeSED(distance=distance)
    dustyMod.PlotSED(xlim=(0,20),ylim=(0,6e4), ax=ax, kwargs=kwargs_dusty)
    plt.show()


    # fit = DustyFit.DustyFit(dustyMod, Data=Dat)

    # Param = {'Temp1':{'theta0':3500,'minimum':2000,'maximum':4000}, 
    #      'Opacity':{'theta0':1,'minimum':0.1,'maximum':4},
    #      'Lest':{'theta0':1,'minimum':0,'maximum':10}, #Mandatory, fit the Lestimation (Luminosity = Lest*Lestimation)
    # }

    # #Initialize the parameter of the MCMC
    # ParamFit = {
    #                         'nsimu': 100,         #Number of iteration (larger is better but slow due to the dusty runtime)   
    #                         'updatesigma': True,  #Update the sigma of the likelihood function
    #                         'method': 'dram',     #Sampling method
    #                         'adaptint': 100,      #Number of interates between adaptation.
    #                         'verbosity': 0,       #0=No output, 1=Output
    #                         'waitbar': True,      #Display a waitbar
    #                     }
    # fit.set_Param(Param)
    # fit.set_ParamFit(ParamFit)

    # fit.LunchFit()
    # fit.get_Fit().PrintResults()
    # fit.PlotResults(xlim=(0,5),ylim=(0,100), ax=ax, kwargs_fit=kwargs_dusty, kwargs_data=kwargs_data)



