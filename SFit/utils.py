import matplotlib.pyplot as plt
import numpy as np
import subprocess
import glob
import pandas as pd
import astropy
import os
import astropy.units as u
from scipy.interpolate import make_interp_spline
from pymcmcstat.MCMC import MCMC



def ScatterPlot(Flux, Wavelength, unit = {'x':'Wavelength', 'y': 'Flux'}, xlim = {} , ylim = {}, scale='linear',kwargs={}):
    """Create a scatter plot

    Args:
        Flux (array): Flux 
        Wavelength (array): Wavelength
        unit (dict, optional): Description of what is show on each axes. Defaults to {'x':'Wavelength', 'y': 'Flux'}.
        xlim (dict, optional): x axis limit.
        ylim (dict, optional): y axis limit.
        scale (str, optional): Y axe scale. Defaults to 'linear'.
        kwargs (dict, optional): matplolib.pyplot kwargs. Defaults to {}.

    """
    fig,ax = plt.subplots()

    ax.scatter(Wavelength,Flux,**kwargs)
    ax.set_xlabel(f'{unit['x']}')
    ax.set_ylabel(f'{unit['y']}')

    ax.set_yscale(scale)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.invert_yaxis()

    plt.show()

def Plot(Flux, Wavelength, unit = {'x':'Wavelength', 'y': 'Flux'}, xlim = {} , ylim = {}, scale='linear',kwargs={}):
    """Create a scatter plot

    Args:
        Flux (array): Flux 
        Wavelength (array): Wavelength
        unit (dict, optional): Description of what is show on each axes. Defaults to {'x':'Wavelength', 'y': 'Flux'}.
        xlim (dict, optional): x axis limit.
        ylim (dict, optional): y axis limit.
        scale (str, optional): Y axe scale. Defaults to 'linear'.
        kwargs (dict, optional): matplolib.pyplot kwargs. Defaults to {}.

    """
    fig,ax = plt.subplots()

    ax.plot(Wavelength,Flux,**kwargs)
    ax.set_xlabel(f'{unit['x']}')
    ax.set_ylabel(f'{unit['y']}')

    ax.set_yscale(scale)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    #ax.invert_yaxis()

    plt.show()

def PrintFile(file, stop=-1):
    """Print a file

    Args:
        file (string): Path to the file
        stop (int, optional): line where to stop printing. Defaults to -1.
    """
    with open(file , 'r') as f:
            lines = f.readlines()
            for line in lines[:stop]:
                print(line)

def LoadFile(file):
    """Load a file in an array

    Args:
        file (string): Path to the file 

    Returns:
        array: array containing the lines of the file
    """

    with open(file , 'r') as f:
        lines = f.readlines()
        return lines
    
def LoadCSV(file, sep=','):
    """Load a csv file

    Args:
        file (string): Path to the file 

    Returns:
        Dataframe: dataframe containing the csv file
    """
    return pd.read_csv(file, sep=sep)

def LoadFits(file):
    """Load a fits file

    Args:
        file (string): Path to the file 

    Returns:
        array: array containing the file
    """
    return astropy.io.fits.open(file)[0].data 

def SearchLine(file,line):
    for i,lines in enumerate(file):
        if line in lines:
            return i 
    raise Exception(f'This line does dot exist: {line}')

# def SearchStr(line,car):
#     for i,seq in enumerate(line.split(' ')):
#         if car in seq:
#             return i
#     return None

# def ModifyStr(line, index, car):
#     l = line.split(' ')
#     l[index] = car
#     return ' '.join(l)

def SaveFile(Path,file):
    with open(Path , 'w') as f:
            f.write("".join(file))

def ChangeParameter(Path,change,car):
     
    file = LoadFile(Path)
    for param in change.keys():
        line = SearchLine(file,car[param])
        new_line = change[param]
        file[line] = new_line

    SaveFile(Path,file)

def ListToDict(keys,values):
     key_value_pairs = zip(keys, values)
     return dict(key_value_pairs)

def SuppCarList(list,car):
     return [el.split('\n')[0] for el in list if el not in car]

def GetColumnSpectrum(file, index, index_header=0):
     array = np.asarray(file[index_header:])
     col = np.array([el.split('  ')[1:-1] for el in array],dtype=float).T[index]
     return col

def WattToJansky(Flux, Wavelength):
     return (Flux * u.W/u.m**2).to(u.Jy,equivalencies=u.spectral_density(Wavelength * u.um))

def UmToMeter(Wavelength):
     return (Wavelength * u.um).to(u.m)

def CalculRayonVrai(results,L):
     return results['r1(cm)']*np.sqrt(L/1e4) *1e-2

def CalculFluxTotal(F,r_vrai,distance):
     return F*(4*np.pi*r_vrai**2)/(4*np.pi*distance**2)

def Interpolate(Wavelength, Flux, order=3):
    return make_interp_spline(Wavelength, Flux, order)

def model(q,x,y):
    I = q
    try:
        return I*Interpolate(x.flatten(),y.flatten())(x.flatten())
    except:
        return I*Interpolate(x,y)(x)

def Chi2(theta,data):
    try:
        xdata = data.xdata[0]
        ydata = data.ydata[0]
    except AttributeError:
        xdata,ydata = data
    
    ymodel = model(theta[:2], xdata, ydata).reshape(ydata.shape)
    # calc sos
    ss = np.nansum((ymodel - ydata)**2)
    return ss

def SetMCMCParam(mc=MCMC, param={}):
     
     for par in param.keys():
          mc.parameters.add_model_parameter(name = par,
                                            theta0 = param[par]['theta0'],
                                            minimum = param[par]['minimum'],
                                            maximum = param[par]['maximum'])
          
def SetMCMC(mc=MCMC,param={}):
    mc.simulation_options.define_simulation_options(nsimu=param['nsimu'],
                                                    updatesigma=param['updatesigma'],
                                                    method=param['method'],
                                                    adaptint=param['adaptint'],
                                                    verbosity=param['verbosity'],
                                                    waitbar=param['waitbar'])

if __name__=="__main__":
     pass
    #  f = LoadFile('/Users/gabriel/Documents/SFit/AFGL4106.inp')
    #  i = SearchLine(f, 'Number of additional components')
    #  j = SearchStr(f[i], '.nk')  
    #  print(f[i])
    #  print(ModifyStr(f[i], j, 'prout'))