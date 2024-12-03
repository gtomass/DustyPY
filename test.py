import DustyPY.stars as stars
import DustyPY.dusty as dusty
import DustyPY.Data as Data
import DustyPY.DustyFit as DustyFit
import DustyPY.utils as utils
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson

if __name__ == '__main__':
    

    S1 = stars.Star('E1',3400,1)
    S2 = stars.Star('E2',7000,1)


    Composition = {'Mg2SiO4':0.04,
                    'Al2O3-comp':0.15,
                    'FeO':0.01,
                    'MgSiO3':0.1,
                    'H2Oice':0.05,
                    'MgFeSiO4':0.65}

    # Composition = {'sil-draine':1}

    
    # density = {'density type':'RDWA', 'shell': 1000}
    density = None

    dust = stars.Dust(DustSize={'Distribution':'MRN'},Composition = Composition, tau = 1, Sublimation = 1500, Properties='common_and_addl_grain_composite', Temperature=130, Density=density)

    mod = stars.Model('AFGL4106', NbStar=2, Stars=[S1, S2],
                        dust=dust, distance=1470,
                        Spectral='black_body',
                        SpectralFile='/Users/gabriel/Documents/These/Recherche/data/flux_4000_25000_1.000_0.000.dat',
                        SiOAbsorption=10)

    dustyMod = dusty.Dusty(PATH='/Users/gabriel/Documents/These/Recherche/lib/dusty',
                            model=mod, Lestimation=1e4)

    

    fig, ax = plt.subplots()
    kwargs_data = {'fmt': 'dr', 'markersize':5,'label': 'Data'}
    kwargs_dusty = {'color': 'k', 'label': 'Dusty','linewidth': 1}

    

    Dat = Data.Data()

    table = Dat.querry_vizier_data(radius = 5, target='AFGL4106')

    ISO = Dat.import_data('/Users/gabriel/Documents/These/Recherche/data/24900158_sws.fit').T
    wave, flux = ISO[0], ISO[1]

    band = utils.get_bandpass_name()
    band = [b for b in band if 'iso' in b]
    filters = [utils.get_bandpass(f) for f in band]

    fluxes = utils.integrate_SED_bandpass(wave, flux, dict(zip(band,band)))
    central = [3e8/(utils.get_central_wavelegnth(f)*1e-10) * 1e-9 for f in filters]

    for c, f, b in zip(central, fluxes, band):
        table.add_row(dict(zip(['sed_freq', 'sed_flux', 'sed_filter'], [c, f, b])))

    dustyMod.make_wavelength()

    table = Dat.restrict_data_vizier(table)

    Dat.set_vizier_data(table)
    Dat.restrict_data(['xdata <= 30'])
    Dat.unred_data(EBV=1.07)
    dustyMod.change_parameter()
    dustyMod.lunch_dusty()
    dustyMod.make_SED(distance=1470)

    dustyMod.plot_SED(xlim=(0.1,110), ylim=(.001,3000),scale='log', ax=ax, kwargs=kwargs_dusty)
    Dat.scatter_data(ax=ax, kwargs=kwargs_data)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.show()

    # fit = DustyFit.DustyFit(dustyMod, data=Dat)

    # Param = {'Temp1':{'theta0':3000,'minimum':2900,'maximum':3500},
    #         'Temp2':{'theta0':7000,'minimum':6000,'maximum':8000},
    #         'Lum1':{'theta0':1,'minimum':0,'maximum':3},
    #         'Lum2':{'theta0':1,'minimum':0,'maximum':3},
    #         'Opacity':{'theta0':1,'minimum':0.1,'maximum':2},
    #         'Lest':{'theta0':3,'minimum':0,'maximum':7}} #Mandatory, fit the Lestimation (Luminosity = Lest*Lestimation)

    # #Initialize the parameter of the MCMC
    # ParamFit = {
    #                         'nsimu': 1000,         #Number of iteration (larger is better but slow due to the dusty runtime)
    #                         'updatesigma': True,  #Update the sigma of the likelihood function
    #                         'method': 'dram',     #Sampling method
    #                         'adaptint': 50,      #Number of interates between adaptation.
    #                         'verbosity': 0,       #0=No output, 1=Output
    #                         'waitbar': True,      #Display a waitbar
    #                     }
    # fit.set_Param(Param)
    # fit.set_ParamFit(ParamFit)

    # fit.lunch_fit(chi2='Chi2_modified')
    # fit.get_Fit().print_results()
    # fit.plot_stats()
    # ax.set_xscale('log')
    # fit.plot_results(xlim=(0.1,110), ylim=(.001,3000),scale='log', ax=ax, kwargs_fit=kwargs_dusty, kwargs_data=kwargs_data, save=True,
    #                  unit={'x': r'Wavelength ($\mu$m)', 'y': 'Flux (Jy)'})

