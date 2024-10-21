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
from PyAstronomy import pyasl



def ScatterPlot(Flux, Wavelength, unit=None, xlim=None, ylim=None, ax = None, scale='linear', kwargs=None):
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
    if unit is None:
        unit = {'x':'Wavelength', 'y': 'Flux'}
    if xlim is None:
        xlim = {}
    if ylim is None:
        ylim = {}
    if kwargs is None:
        kwargs = {}
    if ax is None:
        fig,ax = plt.subplots()

    ax.scatter(Wavelength,Flux,**kwargs)
    ax.set_xlabel(unit['x'])
    ax.set_ylabel(unit['y'])

    ax.set_yscale(scale)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    #ax.invert_yaxis()

    #plt.show()

def Plot(Flux, Wavelength, unit=None, xlim=None, ylim=None, ax=None, scale='linear', kwargs=None):
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
    if unit is None:
        unit = {'x':'Wavelength', 'y': 'Flux'}
    if xlim is None:
        xlim = {}
    if ylim is None:
        ylim = {}
    if kwargs is None:
        kwargs = {}
    if ax is None:
        fig,ax = plt.subplots()

    ax.plot(Wavelength,Flux,**kwargs)
    ax.set_xlabel(f'{unit['x']}')
    ax.set_ylabel(f'{unit['y']}')

    ax.set_yscale(scale)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)


    #plt.show()

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

def LoadFile(Path):
    """Load a file in an array

    Args:
        Path (string): Path to the file 

    Returns:
        array: array containing the lines of the file
    """

    with open(Path , 'r') as f:
        return f.readlines()
    
def LoadCSV(Path, sep=','):
    """Load a csv file

    Args:
        Path (string): Path to the file 

    Returns:
        Dataframe: dataframe containing the csv file
    """
    return pd.read_csv(Path, sep=sep)

def LoadFits(Path):
    """Load a fits file

    Args:
        Path (string): Path to the file 

    Returns:
        array: array containing the file
    """
    return astropy.io.fits.open(Path)[0].data 


def SearchLine(file, line):
    """
    Searches for a specific line in a given file and returns the index of the line if found.

    Args:
        file (iterable): An iterable object representing the file, where each element is a line in the file.
        line (str): The line to search for within the file.

    Returns:
        int: The index of the line if found.

    Raises:
        Exception: If the specified line does not exist in the file.
    """
    for i,lines in enumerate(file):
        if line in lines:
            return i 
    raise Exception(f'This line does dot exist: {line}')

def DeleteLine(file, line):
    """
    Deletes a specific line in a given file.

    Args:
        file (iterable): An iterable object representing the file, where each element is a line in the file.
        line (str): The line to delete within the file.

    Returns:
        iterable: The file with the specified line removed.

    Raises:
        Exception: If the specified line does not exist in the file.
    """
    for i,lines in enumerate(file):
        if line in lines:
            del file[i]
            return file
    raise Exception(f'This line does dot exist: {line}')


def SaveFile(Path,file):
    with open(Path , 'w') as f:
            f.write("".join(file))

def ChangeParameter(Path, change, car, ncomp):
    file = LoadFile(Path)
    
    # Remove all lines containing '.nk'
    file = [line for line in file if ('.nk' not in line) and ('cocomponents' not in line)] 
    
    for param in change.keys():
        line = SearchLine(file, car[param])
        new_line = change[param]
        file[line] = new_line

    SaveFile(Path, file)

def ListToDict(keys,values):
     key_value_pairs = zip(keys, values)
     return dict(key_value_pairs)

def SuppCarList(list,car):
     return [el.split('\n')[0] for el in list if el not in car]

def GetColumnSpectrum(file, index, index_header=0):
    array = np.asarray(file[index_header:])
    return np.array([el.split('  ')[1:-1] for el in array],dtype=float).T[index]

def WattToJansky(Flux, Wavelength):
     #return Flux *1e26*(Wavelength)**2/3e8
     return (Flux * u.W/u.m**2).to(u.Jy,equivalencies=u.spectral_density(Wavelength * u.um)).value

def JanskyToWatt(Flux, Wavelength):
     #return Flux * 1e-26*3e8/Wavelength**2
     return (Flux * u.Jy ).to(u.W/u.m**2,equivalencies=u.spectral_density(Wavelength * u.um)).value

def UmToMeter(Wavelength):
     return (Wavelength * u.um).to(u.m).value

def CalculRayonVrai(results,L):
     return results['r1(cm)']*np.sqrt(L/1e4) *1e-2

def CalculFluxTotal(F,r_vrai,distance):
     return F*(4*np.pi*r_vrai**2)/(4*np.pi*distance**2)

def Interpolate(Wavelength, Flux, order=3):
    return make_interp_spline(Wavelength, Flux, order)

def model(q,xdata,xdusty,ydusty):
    I = q
    try:
        return I*Interpolate(np.asarray(xdusty).flatten(),np.asarray(ydusty).flatten())(np.asarray(xdata).flatten())
    except Exception:
        return I*Interpolate(xdusty,ydusty)(xdusty)

def Chi2(theta,data):
    try:
        xdata = data.xdata[0]
        ydata = data.ydata[0]
    except AttributeError:
        xdata,ydata = data

    ymodel = model(theta, xdata, ydata).reshape(ydata.shape)
    return np.nansum((ymodel - ydata)**2)

def SetMCMCParam(mc=MCMC, param=None):
     
    if param is None:
        param = {}
    for par in param.keys():
         mc.parameters.add_model_parameter(name = par,
                                           theta0 = param[par]['theta0'],
                                           minimum = param[par]['minimum'],
                                           maximum = param[par]['maximum'])
          
def SetMCMC(mc=MCMC,param=None):
    mc.simulation_options.define_simulation_options(nsimu=param['nsimu'],
                                                    updatesigma=param['updatesigma'],
                                                    method=param['method'],
                                                    adaptint=param['adaptint'],
                                                    verbosity=param['verbosity'],
                                                    waitbar=param['waitbar'])


def Unred(Wavelength, Flux, EBV, Rv=3.1):
    return pyasl.unred(Wavelength,Flux, ebv=EBV, R_V=Rv)

def Check(change):
        L = 0
        for key in change.keys():
            if 'Lum' in key:
                L+=change[key]
        return L==1



if __name__=="__main__":
     pass