from DustyPY import utils as utils
from DustyPY import stars as stars
from DustyPY import dusty as dusty
import matplotlib.pyplot as plt


if __name__=='__main__':
    
    S = stars.Star(Name='E1',Temperature=3500,Luminosity=1)

    dust = stars.Dust(DustSize={'amin' : 0.01, 'amax': 1}, 
                      Composition={'Al2O3-comp':0.2,'MgFeSiO4':0.16,'Ca2Al2SiO7':0.64},
                      tau=1,
                      Sublimation=1200,
                      Properties='common_and_addl_grain_composite'
                      )
    
    mod = stars.Model('Prout', NbStar=1, Stars=[S], 
                      Dust=dust, distance=197.0,
                      Spectral='black_body',
                      SiOAbsorption=0.5)

    
    dustyMod = dusty.Dusty(PATH='/Users/gabriel/Documents/These/Recherche/lib/dusty/release/dusty',
                           Model=mod,Lestimation=1.8e5)
    
    print(dustyMod.AvailableComposition())
    dustyMod.ChangeParameter()
    dustyMod.LunchDusty()
    dustyMod.MakeSED(distance=197)
    fig,ax = plt.subplots()
    dustyMod.PlotSED(xlim=(0,20),ax=ax)
    plt.show()
