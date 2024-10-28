from DustyPY import utils as utils
from DustyPY import stars as stars
from DustyPY import dusty as dusty
from DustyPY import DustyFit as DustyFit
from DustyPY import Data as Data
import matplotlib.pyplot as plt


if __name__=='__main__':
    
    S = stars.Star(Name='E1',Temperature=3697.87,Luminosity=1)

    dust = stars.Dust(DustSize={'Distribution':'MODIFIED_MRN','amin' : 0.01, 'amax': 1}, 
                      Composition={'Al2O3-comp':0.2,'MgFeSiO4':0.16,'Ca2Al2SiO7':0.64},
                      tau=1.06,
                      Sublimation=1200,
                      Properties='common_and_addl_grain_composite'
                      )
    
    mod = stars.Model('Betelgeuse', NbStar=1, Stars=[S], 
                      Dust=dust, distance=197.0,
                      Spectral='black_body',
                      SiOAbsorption=10)

    
    dustyMod = dusty.Dusty(PATH='/Users/gabriel/Documents/These/Recherche/lib/dusty/release/dusty',
                           Model=mod,Lestimation=1.8e5)
    
    dustyMod.MakeWavelength(119)

    fig,ax = plt.subplots()
    kwargs_data = {'fmt':'+','color':'r'}
    kwargs_dusty = {'color':'k'}

    Dat = Data.Data()

    table = Dat.QuerryVizierData(radius = 5, target='Betelgeuse')
    Dat.SetVizierData(table)
    Dat.RestrictData(['yerr != 0.'])

    Dat.UnredData(EBV=0.4042)
    
    # dustyMod.ChangeParameter()
    # dustyMod.LunchDusty(verbose=0)
    # dustyMod.MakeSED(distance=197)
    # dustyMod.PlotSED(xlim=(0,10),ax=ax,kwargs=kwargs_dusty)
    # Dat.ScatterData(xlim=(0,20), ax=ax, kwargs=kwargs_data)
    # plt.show()

    

    fit = DustyFit.DustyFit(dustyMod, Data=Dat)

    Param = {'Temp1':{'theta0':3500,'minimum':2000,'maximum':4000}, 
         'Opacity':{'theta0':1,'minimum':0.1,'maximum':4},
         'Lest':{'theta0':1,'minimum':0,'maximum':10}, #Mandatory, fit the Lestimation (Luminosity = Lest*Lestimation)
    }

    #Initialize the parameter of the MCMC
    ParamFit = {
                            'nsimu': 1000,         #Number of iteration (larger is better but slow due to the dusty runtime)   
                            'updatesigma': True,  #Update the sigma of the likelihood function
                            'method': 'dram',     #Sampling method
                            'adaptint': 100,      #Number of interates between adaptation.
                            'verbosity': 0,       #0=No output, 1=Output
                            'waitbar': True,      #Display a waitbar
                        }
    fit.set_Param(Param)
    fit.set_ParamFit(ParamFit)

    fit.LunchFit()
    fit.get_Fit().PrintResults()
    fit.PlotResults(xlim=(0,5), ax=ax, kwargs_fit=kwargs_dusty, kwargs_data=kwargs_data)

