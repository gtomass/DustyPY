import SFit.stars as stars
import SFit.dusty as dusty
import SFit.utils as utils

if __name__== '__main__':

    distance = 1000*3.26*365*86400*3e8

    #print(7.6/(7.6+9.8))

    S1 = stars.Star(Name='E1',Temperature=3000,Luminosity=0.57)
    S2 = stars.Star(Name='E2',Temperature=7350,Luminosity=0.43)

    Composition = {'Mg2SiO4':0.04, 
                   'MgSiO3':0.1, 
                   'H2Oice': 0.05,
                   'Al2O3-comp': 0.15,
                   'MgFeSiO4':0.65,
                   'FeO':0.01 
                   }

    dust = stars.Dust(DustSize={'amin' : 0.3, 'amax': 10}, Composition=Composition)

    mod = stars.Model('AFGL4106', NbStar=2, Stars=[S1,S2], Dust=dust)
   
    fit = dusty.Dusty(PATH='/Users/gabriel/Documents/Stage/code/dustyV2/', Model=mod)

    #fit.ChangeParameter()
    #fit.LunchDusty()
    print(fit.GetResults())

    kwargs = {'marker':'+', 'color':'red'}

    W,F = fit.GetSED(distance = distance, Jansky=False)

    utils.ScatterPlot(F,W,xlim={0.1,10},ylim={100,0.001},kwargs=kwargs)