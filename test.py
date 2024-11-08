import DustyPY.stars as stars
import DustyPY.dusty as dusty
import DustyPY.Data as Data
import DustyPY.DustyFit as DustyFit
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson

if __name__ == '__main__':

    space = 3

    data = np.loadtxt('/Users/gabriel/Documents/These/Recherche/data/flux_4000_25000_1.000_0.000.dat', skiprows=2,dtype=str).T
    data_copie = np.zeros(data.shape)

    for i in range(len(data)):
        line = data[i]
        for j in range(len(line)):
            data_copie[i][j] = float('E'.join(data[i][j].split('D'))) if 'D' in line[j] else float(line[j])

    nL = len(data_copie[0][::space])
    print(nL)
    F = data_copie[1][::space]

    Fbol = simpson(y = data_copie[1][::space], x = data_copie[0][::space])
    

    S1 = stars.Star('E1',3400,1)


    Composition = {'Mg2SiO4':0.04,
                    'Al2O3-comp':0.15,
                    'FeO':0.01,
                    'MgSiO3':0.1,
                    'H2Oice':0.05,
                    'MgFeSiO4':0.65}

    dust = stars.Dust(Composition = Composition, tau = 1, Sublimation = 1500, Properties='common_and_addl_grain_composite', Temperature=800)

    mod = stars.Model('test', NbStar=1, Stars=[S1],
                        Dust=dust, distance=1000,
                        Spectral='FILE_F_LAMBDA',
                        SpectralFile='/Users/gabriel/Documents/These/Recherche/data/flux_4000_25000_1.000_0.000.dat',
                        SiOAbsorption=10)

    dustyMod = dusty.Dusty(PATH='/Users/gabriel/Documents/These/Recherche/lib/dusty',
                            Model=mod, Lestimation=1e3/3)

    dustyMod.make_wavelength(intervals=[(1e-2,np.min(data_copie[0]), 50), (np.min(data_copie[0]),np.max(data_copie[0]), nL), (np.max(data_copie[0]),1e2, 100)])

    fig, ax = plt.subplots()
    kwargs_data = {'color': 'r', 'label': 'input', 'linewidth': 0.5}
    kwargs_dusty = {'color': 'k', 'label': 'output','linewidth': 0.5}

    

    # Dat = Data.Data()


    # table = Dat.querry_vizier_data(radius = 5, target='AFGL4106')
    # Dat.set_vizier_data(table)
    # Dat.restrict_data(['yerr != 0.'])

    # Dat.unred_data(EBV=1.07)

    dustyMod.change_parameter()
    dustyMod.lunch_dusty(verbose=0)

    ax.plot(data_copie[0][::space], F/np.max(F) , **kwargs_data)
    dustyMod.make_SED(distance=197, normalize=True)
    dustyMod.plot_SED(xlim = (0.2,5),ylim = (-0.1,1.1),ax=ax, kwargs=kwargs_dusty, normalize=True, unit={'x': r'Wavelength ($\mu$m)', 'y': 'Normalized Flux'})
    
   
    # Dat.scatter_data(xlim=(0,10), ax=ax, kwargs=kwargs_data)
    ax.legend()
    #plt.savefig('/Users/gabriel/Documents/These/Recherche/Plot/CO5BOLT_dusty.png', dpi=300)
    plt.show()

    
    # fit = DustyFit.DustyFit(dustyMod, data=Dat)

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

    # fit.lunch_fit()
    # fit.get_Fit().print_results()
    # fit.plot_results(xlim=(0,5), ax=ax, kwargs_fit=kwargs_dusty, kwargs_data=kwargs_data)
