import DustyPY.stars as stars
import DustyPY.dusty as dusty
import DustyPY.Data as Data
import DustyPY.DustyFit as DustyFit
import DustyPY.utils as utils
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson
import astropy.units as u

if __name__ == '__main__':
    def as_si(x, ndp):
        s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
        m, e = s.split('e')
        return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))


    #----------------
    L1 = 7.95
    L2 = 3.85
    T1 = 5955.56
    T2 = 2596.52
    T_dust = 158 #±7.53

    Ldusty = 1e4
    Lest = L1+L2
    Ltot = Ldusty*Lest
    tau = 2. #±0.53

    amin = 0.005
    amax = 0.25
    q = 3.5
    EBV = 1.07

    S1 = stars.Star('E1',T1,L1)
    S2 = stars.Star('E2',T2,L2)

    Jansky = False

    kwargs_data = {'fmt': 'dr', 'markersize':4,'label': 'Data'}
    kwargs_dusty = {'color': 'k', 'label': 'Dusty','linewidth': 1}
    unit = {'x':r'$\lambda$ [$\mu m$]','y':r'$F_\lambda$ [Jy]'} if Jansky else {'x':r'$\lambda$ [$\mu m$]','y':r'$\lambda F_\lambda$ [Wm$^{-2}$]'}



    xlim = (0.1, 170)
    ylim = (1e-13, 1e-8) if not Jansky else (0.1,3000)

    textstr = '\n'.join((
        r'T$_{1}=%.0f$ K' % (T1, ),
        r'T$_{2}=%.0f$ K' % (T2, ),
        r'L$_{\mathrm{tot}}=%s$ $L_{\odot}$' % (as_si(Lest*Ldusty, 2), ),
        r'$T_{\mathrm{dust}}=%.0f$ K' % (T_dust, ),
        r'$\tau=%.2f$' % (tau, ),
        ))

    props = dict(boxstyle='round', facecolor='grey', alpha=0)

    #-----------------

    Composition = {'Mg2SiO4':0.04,
                'MgSiO3':0.1,
                'Al2O3-comp':0.15,
                'MgFeSiO4':0.65,
                'H2Oice':0.05,
                'FeO':0.01}
    density = {'density type':'RDW', 'shell': 40}
    Dustsize = {'Distribution':'MODIFIED_MRN', 'amin':amin, 'amax':amax, 'q':q}
    dust = stars.Dust(DustSize=Dustsize,Composition = Composition, tau = tau, Sublimation = 1500, Properties='common_and_addl_grain_composite', Temperature=T_dust, Density=density)


    mod = stars.Model('AFGL4106', NbStar=2, Stars=[S1,S2],
                            dust=dust, distance=3200,
                            Spectral='black_body')

    dustyMod = dusty.Dusty(PATH='/Users/gtomassini/NextCloud/These/Recherche/lib/dusty',
                                model=mod, Lestimation=Ldusty)

    Dat = Data()

    table = Dat.querry_vizier_data(radius = 5, target='AFGL4106')

    ISO = Dat.import_data('/Users/gtomassini/NextCloud/These/Recherche/data/24900158_sws.fit').T
    wave, flux = ISO[0], ISO[1]

    band = utils.get_bandpass_name()
    band = [b for b in band if 'iso' in b]
    filters = [utils.get_bandpass(f) for f in band]

    fluxes,fluxes_err = utils.integrate_SED_bandpass(wave, flux, dict(zip(band,band)))
    central = [3e8/(utils.get_central_wavelegnth(f)*1e-10) * 1e-9 for f in filters]

    for c, f,fe, b in zip(central, fluxes, fluxes_err, band):
        table.add_row(dict(zip(['_tabname','sed_freq', 'sed_flux', 'sed_eflux', 'sed_filter'], ['ISO SWS',c, f,fe ,b])))

    table = Dat.restrict_data_vizier(table)
    table.remove_rows(np.argwhere(table['sed_filter'] == 'GAIA/GAIA3:Gbp').flatten())
    table.remove_rows(np.argwhere(table['_tabname'] == 'I/353/gsc242').flatten())
    table.remove_rows(np.argwhere(table['_tabname'] == 'II/328/allwise').flatten())
    table.remove_rows(np.argwhere(table['_tabname'] == 'II/366/catv2021').flatten())
    table.remove_rows(np.argwhere(table['_tabname'] == 'IV/38/tic').flatten())
    table.remove_rows(np.argwhere(table['_tabname'] == 'IV/39/tic82').flatten())


    Dat.set_vizier_data(table)
    # Dat.restrict_data(['xdata <= 100'])
    Dat.unred_data(EBV=EBV)
    Dat.convert_to_watt()

    if False:
        fig, ax = plt.subplots()
        Dat.scatter_data(ax=ax,xlim=xlim, ylim=ylim, kwargs=kwargs_data, unit=unit)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid()
        plt.savefig('/data/home/gtomassini/gtomassini/SEDfit/plot/AFGL4106_data.png', dpi=300, bbox_inches='tight')

    dustyMod.make_wavelength(intervals=[(0.01, .4,20),(.401,.8,40),(0.801,20,100),(20.1,100,70),(101,1000,10)])
    # dustyMod.change_parameter()
    # dustyMod.lunch_dusty(verbose=0, logfile=True)


    if False:

        fig, ax = plt.subplots()

        comon_filter = Dat.get_common_filters()

        kwargs_data = {'fmt': 'dr', 'markersize':5,'label': 'Data'}

        dustyMod.make_SED(luminosity=Lest, Jansky=Jansky)
        dustyMod.plot_SED(xlim=xlim, ylim=ylim,scale='log', ax=ax, kwargs=kwargs_dusty, unit=unit, normalize=False)
        dustyMod.get_SED().scatter_SED_bandpass(ax=ax, kwargs={'color':'b','label':'Fit'},scale='log', unit=unit, bandpass = dict(zip(comon_filter.values(),comon_filter.values())))
        Dat.scatter_data(ax=ax,xlim=xlim, ylim=ylim, kwargs=kwargs_data, unit=unit, normalize=False)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid()
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=11,
                verticalalignment='top', bbox=props)
        plt.savefig('/data/home/gtomassini/gtomassini/SEDfit/plot/AFGL4106_dusty.png', dpi=300)


    fit = DustyFit(dustyMod, data=Dat, ncpu=4)

    Param = {
                'Temp1':{'theta0':T1,'minimum':T1-1000,'maximum':T1+1000, 'sample':False},
                'Temp2':{'theta0':T2,'minimum':T2-1000,'maximum':T2+1000, 'sample':False},
                'Lum1':{'theta0':L1,'minimum':L1-0.2,'maximum':L1+0.2, 'sample':True},
                'Lum2':{'theta0':L2,'minimum':L2-0.2,'maximum':L2+0.2, 'sample':True},
                'Opacity':{'theta0':tau,'minimum':1,'maximum':5, 'sample':True},
                'Temperature':{'theta0':T_dust,'minimum':120,'maximum':180, 'sample':True},
                'amin':{'theta0':np.log10(amin),'minimum':np.log10(0.005),'maximum':np.log10(0.4), 'sample':False},
                'amax':{'theta0':np.log10(amax),'minimum':np.log10(0.15),'maximum':np.log10(6), 'sample':False},
                'q':{'theta0':q,'minimum':2.5,'maximum':3.5, 'sample':False},
                'Lest':{'theta0':Lest,'minimum':Lest-0.4,'maximum':Lest+0.4, 'sample':True},
                'shell':{'theta0':40,'minimum':20,'maximum':60, 'sample':True},
                } 
    fit.set_Param(Param)

    ParamFit = {
                            'nsimu': 250,         #Number of iteration (larger is better but slow due to the dusty runtime)
                            'updatesigma': True,  #Update the sigma of the likelihood function
                            'method': 'dram',     #Sampling method
                            'adaptint': 10,      #Number of interates between adaptation.
                            'verbosity': 0,       #0=No output, 1=Output
                            'waitbar': True,      #Display a waitbar
                        }
    fit.set_ParamFit(ParamFit)

    fit.lunch_fit(chi2='Chi2_modified', logfile=True, Jansky=Jansky)

    stats = fit.get_Fit().get_Stats()['mean']

    p = [key for key, value in Param.items() if value['sample']]
    T1_fit = stats[p.index('Temp1')] if Param['Temp1']['sample'] else Param['Temp1']['theta0']
    T2_fit = stats[p.index('Temp2')] if Param['Temp2']['sample'] else Param['Temp2']['theta0']
    Lest_fit = stats[p.index('Lest')] if Param['Lest']['sample'] else Param['Lest']['theta0']
    L1_fit = stats[p.index('Lum1')] if Param['Lum1']['sample'] else Param['Lum1']['theta0']
    L2_fit = stats[p.index('Lum2')] if Param['Lum2']['sample'] else Param['Lum2']['theta0']
    T_dust_fit = stats[p.index('Temperature')] if Param['Temperature']['sample'] else Param['Temperature']['theta0']
    tau_fit = stats[p.index('Opacity')] if Param['Opacity']['sample'] else Param['Opacity']['theta0']
    amin_fit = np.round(10**stats[p.index('amin')],3) if Param['amin']['sample'] else np.round(10**Param['amin']['theta0'],3)
    amax_fit = np.round(10**stats[p.index('amax')],3) if Param['amax']['sample'] else np.round(10**Param['amax']['theta0'],3)
    q_fit = stats[p.index('q')] if Param['q']['sample'] else Param['q']['theta0']
    shell_fit = stats[p.index('shell')] if Param['shell']['sample'] else Param['shell']['theta0']

    chi2 = fit.get_chi2(chi2='Chi2_modified', Jansky=Jansky)

    textstr_fit = '\n'.join((
        r'$\chi^{2}$ = %.3f' % (chi2, ),
        r'T$_{1}=%.0f$ K' % (T1_fit, ),
        'T$_{2}=%.0f$ K' % (T2_fit, ),
        r'L$_{\mathrm{tot}}=%s$ $L_{\odot}$' % (as_si(Lest_fit*Ldusty, 2), ),
        r'$T_{\mathrm{dust}}=%.0f$ K' % (T_dust_fit, ),
        r'$\tau=%.2f$' % (tau_fit, ),
        ))

    fig,ax = plt.subplots()
    ax.set_xscale('log')
    ax.grid()
    ax.text(0.03, 0.97, textstr_fit, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    fit.plot_results(xlim=xlim, ylim=ylim,scale='log', ax=ax, kwargs_fit=kwargs_dusty, kwargs_data=kwargs_data, save=True,
                    unit=unit, Jansky=Jansky)
    plt.savefig('/data/home/gtomassini/gtomassini/SEDfit/plot/AFGL4106_fit.png', dpi=300, bbox_inches='tight')
