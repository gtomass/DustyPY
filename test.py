import DustyPY.stars as stars
import DustyPY.dusty as dusty
import DustyPY.Data as Data
import DustyPY.DustyFit as DustyFit
import DustyPY.utils as utils
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson

if __name__ == '__main__':
    L1 = 0.73
    L2 = 0.65
    T1 = 4555
    T2 = 6630
    T_dust = 151.78 #±7.53

    Ldusty = 1e4
    Lest = 2.17
    Ltot = Ldusty*Lest
    tau = 2.01 #±0.53

    S1 = stars.Star('E1',3400,1)
    S2 = stars.Star('E2',7000,1)


    Composition = {'Mg2SiO4':0.04,
                    'Al2O3-comp':0.15,
                    'FeO':0.01,
                    'MgSiO3':0.1,
                    'H2Oice':0.05,
                    'MgFeSiO4':0.65}

    # Composition = {'sil-draine':1}

    
    density = {'density type':'RDWA', 'shell': 1000}
    #density = None

    dust = stars.Dust(DustSize={'Distribution':'MODIFIED_MRN'},Composition = Composition, tau = 1, Sublimation = 1500, Properties='common_and_addl_grain_composite', Temperature=130, Density=density)

    mod = stars.Model('AFGL4106', NbStar=2, Stars=[S1, S2],
                        dust=dust, distance=1470,
                        Spectral='black_body',
                        SpectralFile='/Users/gabriel/Documents/These/Recherche/data/flux_4000_25000_1.000_0.000.dat',
                        SiOAbsorption=10)

    dustyMod = dusty.Dusty(PATH='/Users/gtomassini/NextCloud/These/Recherche/lib/dusty',
                            model=mod, Lestimation=1e4)

    dustyMod.make_wavelength(intervals=[(0.01,0.8,100),(0.801,20,70),(20.01,100,50),(101,1000,10)])
    

    fig, ax = plt.subplots()
    kwargs_data = {'fmt': 'dr', 'markersize':5,'label': 'Data'}
    kwargs_dusty = {'color': 'k', 'label': 'Dusty','linewidth': 1}

    

    Dat = Data.Data()

    table = Dat.querry_vizier_data(radius = 5, target='AFGL4106')

    ISO = Dat.import_data('/Users/gtomassini/NextCloud/These/Recherche/data/24900158_sws.fit').T
    wave, flux = ISO[0], ISO[1]

    band = utils.get_bandpass_name()
    band = [b for b in band if 'iso' in b]
    filters = [utils.get_bandpass(f) for f in band]

    fluxes = utils.integrate_SED_bandpass(wave, flux, dict(zip(band,band)))
    central = [3e8/(utils.get_central_wavelegnth(f)*1e-10) * 1e-9 for f in filters]

    for c, f, b in zip(central, fluxes, band):
        table.add_row(dict(zip(['_tabname','sed_freq', 'sed_flux', 'sed_filter'], ['ISO SWS',c, f, b])))

    table = Dat.restrict_data_vizier(table)
    table.remove_rows(np.argwhere(table['sed_filter'] == 'GAIA/GAIA3:G').flatten())
    table.remove_rows(np.argwhere(table['sed_filter'] == 'WISE:W4').flatten())
    table.remove_rows(np.argwhere(table['sed_filter'] == 'GAIA/GAIA3:Grp').flatten())
    Dat.set_vizier_data(table)
    # Dat.restrict_data(['xdata <= 100'])
    Dat.unred_data(EBV=1.07)

    # Dat.write_table_to_latex(Path = '/Users/gabriel/Documents/These/Recherche/data/AFGL4106.txt', columns=['Wavelength', 'sed_flux', 'sed_eflux', '_tabname'], 
    #                          column_names=['Wavelength', 'Flux', 'Error', 'Catalog'], wavelength=True)
    # dustyMod.change_parameter()
    # dustyMod.lunch_dusty()
    # dustyMod.make_SED(distance=1470)

    # dustyMod.plot_SED(xlim=(0.1,110), ylim=(.001,3000),scale='log', ax=ax, kwargs=kwargs_dusty)
    # Dat.scatter_data(ax=ax, kwargs=kwargs_data)
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # ax.set_xlim(0.1,110)
    # ax.set_ylim(.1,3000)
    # plt.show()

    fit = DustyFit.DustyFit(dustyMod, data=Dat)

    
    Param = {
            'Temp1':{'theta0':T1,'minimum':3800,'maximum':5800, 'sample':True},
            'Temp2':{'theta0':T2,'minimum':5800,'maximum':7800, 'sample':True},
            'Lum1':{'theta0':L1,'minimum':0.1,'maximum':1, 'sample':True},
            'Lum2':{'theta0':L2,'minimum':0.1,'maximum':1, 'sample':True},
            'Opacity':{'theta0':tau,'minimum':1,'maximum':3, 'sample':True},
            'Temperature':{'theta0':T_dust,'minimum':120,'maximum':180, 'sample':True},
            'shell':{'theta0':50,'minimum': 10,'maximum':100, 'sample':True},
            'Lest':{'theta0':Lest,'minimum':.5,'maximum':5, 'sample':True},
            'amin':{'theta0':np.log10(0.005),'minimum':np.log10(0.005),'maximum':np.log10(0.4), 'sample':False},
            'amax':{'theta0':np.log10(0.25),'minimum':np.log10(0.1),'maximum':np.log10(1), 'sample':False},
            'q':{'theta0':3.5,'minimum':2.8,'maximum':3.5, 'sample':False}
            
            } 

    #Initialize the parameter of the MCMC
    ParamFit = {
                            'nsimu': 100,         #Number of iteration (larger is better but slow due to the dusty runtime)
                            'updatesigma': True,  #Update the sigma of the likelihood function
                            'method': 'dram',     #Sampling method
                            'adaptint': 10,      #Number of interates between adaptation.
                            'verbosity': 0,       #0=No output, 1=Output
                            'waitbar': True,      #Display a waitbar
                        }
    fit.set_Param(Param)
    fit.set_ParamFit(ParamFit)

    fit.lunch_fit(chi2='Chi2_modified', logfile=True)
    fit.get_Fit().print_results()
    print(fit.get_chi2(chi2='Chi2_modified'))
    fit.plot_stats()
    fit.plot_results(xlim=(0.1,110), ylim=(.001,3000), scale='log-log', ax=ax, kwargs_fit=kwargs_dusty, kwargs_data=kwargs_data, save=False,
                     unit={'x': r'Wavelength ($\mu$m)', 'y': 'Flux (Jy)'})


    fig = {'xscale': 'log',
            'yscale': 'log', 
            'xlim': (0.1,110), 
            'ylim': (.001,3000), 
            'xlabel': r'Wavelength ($\mu$m)', 
            'ylabel': 'Flux (Jy)', 
            'title': None}
    
    fit.plot_interval(wavelength_intervals=[(0.1,0.8,100),(0.801,20,100),(20.01,100,100),(100.01,1000,10)],
                      ciset={'limits': [50,90,95,99]}, fig=fig)
    plt.show()
    

