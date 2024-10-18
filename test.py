import SFit.stars as stars
import SFit.dusty as dusty

if __name__== '__main__':
    S1 = stars.Star(Name='E1',Temperature=4000,Luminosity=0.5)
    S2 = stars.Star(Name='E2',Temperature=7000,Luminosity=0.5)

    dust = stars.Dust(DustSize={'amin' : 0.3, 'amax': 10}, Composition={'Mg2SiO4':0.8, 'FeO':0.2, "Fe": 0.1})

    mod = stars.Model('AFGL4106', NbStar=2, Stars=[S1,S2], Dust=dust)
   
    fit = dusty.Dusty(PATH='/Users/gabriel/Documents/Stage/code/dustyV2/', Model=mod)

    fit.ChangeParameter()
    fit.LunchDusty()

    # print(fit.get_Model().get_Dust().get_Composition())

    # fit.ChangeParameter()
    # fit.PrintParam()

    #fit.test()


    #fit = Dusty(PATH = '/Users/gabriel/Documents/Stage/code/dustyV2/')

    #fit.ChangeParameter(Star=S,tau=0.2)
    #fit.PrintParam(Star=S)
    #fit.LunchDusty()