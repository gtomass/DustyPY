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

def ErrorPlot(Flux, Wavelength, eFlux, unit=None, xlim=None, ylim=None, ax=None, scale='linear', kwargs=None):
    """Create a scatter plot
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

    ax.errorbar(Wavelength,Flux,yerr=eFlux,**kwargs)
    ax.set_xlabel(f'{unit['x']}')
    ax.set_ylabel(f'{unit['y']}')

    ax.set_yscale(scale)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)


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
    def SaveFile(Path, file):
        """
        Save the given content to a file at the specified path.

        Args:
            Path (str): The path where the file will be saved.
            file (list of str): The content to be written to the file, provided as a list of strings.

        Returns:
            None
        """
    with open(Path , 'w') as f:
            f.write("".join(file))

def ChangeParameter(Path, change, car, ncomp):
    """
    Modify specific parameters in a file by removing certain lines and updating others.
    Args:
        Path (str): The path to the file to be modified.
        change (dict): A dictionary where keys are parameter names and values are the new lines to replace the old ones.
        car (dict): A dictionary where keys are parameter names and values are the search criteria for locating the lines to be changed.
        ncomp (int): The number of components (not used in the current implementation).
    Returns:
        None
    """
    file = LoadFile(Path)
    
    # Remove all lines containing '.nk'
    file = [line for line in file if ('.nk' not in line)] 
    
    for param in change.keys():
        line = SearchLine(file, car[param])
        new_line = change[param]
        file[line] = new_line

    SaveFile(Path, file)

def ListToDict(keys,values):
    """
    Convert two lists into a dictionary.

    This function takes two lists, one containing keys and the other containing values,
    and combines them into a dictionary where each key from the first list is paired
    with the corresponding value from the second list.

    Parameters:
    keys (list): A list of keys.
    values (list): A list of values.

    Returns:
    dict: A dictionary with keys from the 'keys' list and values from the 'values' list.

    Example:
    >>> keys = ['a', 'b', 'c']
    >>> values = [1, 2, 3]
    >>> ListToDict(keys, values)
    {'a': 1, 'b': 2, 'c': 3}
    """
    key_value_pairs = zip(keys, values)
    return dict(key_value_pairs)

def SuppCarList(list,car):
    """
    Filters out elements from the input list that are present in the car list and removes newline characters.

    Args:
        list (list of str): The input list of strings to be filtered.
        car (list of str): The list of strings to be excluded from the input list.

    Returns:
        list of str: A new list with elements from the input list that are not in the car list, 
                        with newline characters removed.
    """
    return [el.split('\n')[0] for el in list if el not in car]

def GetColumnSpectrum(file, index, index_header=0):
    """
    Extrait une colonne spécifique d'un fichier.

    Paramètres:
    file (list of str): Le contenu du fichier en entrée sous forme de liste de chaînes, où chaque chaîne représente une ligne du fichier.
    index (int): L'index de la colonne à extraire.
    index_header (int, optional): L'index de la ligne d'en-tête à ignorer. Par défaut à 0.

    Retourne:
    numpy.ndarray: Un tableau numpy contenant les valeurs de la colonne spécifiée.
    """
    array = np.asarray(file[index_header:])
    return np.array([el.split('  ')[1:-1] for el in array], dtype=float).T[index]

def WattToJansky(Flux, Wavelength):
    """
    Convertit le flux de Watts par mètre carré en Jansky.

    Paramètres:
    Flux (float): Le flux en Watts par mètre carré.
    Wavelength (float): La longueur d'onde en micromètres.

    Retourne:
    float: Le flux en Jansky.
    """
    return (Flux * u.W/u.m**2).to(u.Jy, equivalencies=u.spectral_density(Wavelength * u.um)).value

def JanskyToWatt(Flux, Wavelength):
    """
    Convertit le flux de Jansky en Watts par mètre carré.

    Paramètres:
    Flux (float): Le flux en Jansky.
    Wavelength (float): La longueur d'onde en micromètres.

    Retourne:
    float: Le flux en Watts par mètre carré.
    """
    return (Flux * u.Jy).to(u.W/u.m**2, equivalencies=u.spectral_density(Wavelength * u.um)).value

def UmToMeter(Wavelength):
    """
    Convertit la longueur d'onde de micromètres en mètres.

    Paramètres:
    Wavelength (float): La longueur d'onde en micromètres.

    Retourne:
    float: La longueur d'onde en mètres.
    """
    return (Wavelength * u.um).to(u.m).value

def CalculRayonVrai(results, L):
    """
    Calcule le rayon réel en fonction des résultats et de la luminosité.

    Paramètres:
    results (dict): Un dictionnaire contenant les résultats, y compris 'r1(cm)'.
    L (float): La luminosité en unités arbitraires.

    Retourne:
    float: Le rayon réel en mètres.
    """
    return results['r1(cm)'] * np.sqrt(L / 1e4) * 1e-2

def CalculFluxTotal(F, r_vrai, distance):
    """
    Calcule le flux total en fonction du flux, du rayon réel et de la distance.

    Paramètres:
    F (float): Le flux en unités arbitraires.
    r_vrai (float): Le rayon réel en mètres.
    distance (float): La distance en mètres.

    Retourne:
    float: Le flux total en unités arbitraires.
    """
    return F * (4 * np.pi * r_vrai**2) / (4 * np.pi * distance**2)

def Interpolate(Wavelength, Flux, order=3):
    """
    Interpole les données de flux en fonction de la longueur d'onde.

    Paramètres:
    Wavelength (array-like): Les longueurs d'onde.
    Flux (array-like): Les flux correspondants.
    order (int, optional): L'ordre de l'interpolation. Par défaut à 3.

    Retourne:
    function: Une fonction d'interpolation spline.
    """
    return make_interp_spline(Wavelength, Flux, order)

def model(q, xdata, xdusty, ydusty):
    """
    Calcule le modèle en fonction des paramètres fournis.

    Paramètres:
    q (float): Le facteur d'échelle.
    xdata (array-like): Les données x.
    xdusty (array-like): Les longueurs d'onde du modèle de poussière.
    ydusty (array-like): Les flux du modèle de poussière.

    Retourne:
    array-like: Les valeurs du modèle interpolé.
    """
    I = q
    try:
        return I * Interpolate(np.asarray(xdusty).flatten(), np.asarray(ydusty).flatten())(np.asarray(xdata).flatten())
    except Exception:
        return I * Interpolate(xdusty, ydusty)(xdusty)

def Chi2(theta, data):
    """
    Calcule le chi2 entre les données observées et le modèle.

    Paramètres:
    theta (array-like): Les paramètres du modèle.
    data (tuple or object): Les données observées. Peut être un tuple (xdata, ydata) ou un objet avec des attributs xdata et ydata.

    Retourne:
    float: La valeur du chi2.
    """
    try:
        xdata = data.xdata[0]
        ydata = data.ydata[0]
    except AttributeError:
        xdata, ydata = data

    ymodel = model(theta, xdata, ydata).reshape(ydata.shape)
    return np.nansum((ymodel - ydata)**2)

def SetMCMCParam(mc=MCMC, param=None):
    """
    Définit les paramètres du modèle pour l'objet MCMC.

    Paramètres:
    mc (MCMC object): L'objet MCMC.
    param (dict, optional): Un dictionnaire contenant les paramètres du modèle. Par défaut à None.
    """
    if param is None:
        param = {}
    for par in param.keys():
        mc.parameters.add_model_parameter(name=par,
                                          theta0=param[par]['theta0'],
                                          minimum=param[par]['minimum'],
                                          maximum=param[par]['maximum'])

def SetMCMC(mc=MCMC, param=None):
    """
    Définit les options de simulation pour l'objet MCMC.

    Paramètres:
    mc (MCMC object): L'objet MCMC.
    param (dict, optional): Un dictionnaire contenant les options de simulation. Par défaut à None.
    """
    mc.simulation_options.define_simulation_options(nsimu=param['nsimu'],
                                                    updatesigma=param['updatesigma'],
                                                    method=param['method'],
                                                    adaptint=param['adaptint'],
                                                    verbosity=param['verbosity'])


def Unred(Wavelength, Flux, EBV, Rv=3.1):
    """
    Applique une correction de déreddening aux flux observés.

    Paramètres:
    Wavelength (array-like): Les longueurs d'onde des observations.
    Flux (array-like): Les flux observés.
    EBV (float): La valeur de E(B-V) pour la correction.
    Rv (float, optional): Le rapport de la loi d'extinction. Par défaut à 3.1.

    Retourne:
    array-like: Les flux corrigés.
    """
    return pyasl.unred(Wavelength, Flux, ebv=EBV, R_V=Rv)

def Check(change):
    """
    Vérifie si la somme des changements de luminosité est égale à 1.

    Paramètres:
    change (dict): Un dictionnaire contenant les changements de paramètres.

    Retourne:
    bool: True si la somme des changements de luminosité est égale à 1, sinon False.
    """
    L = 0
    for key in change.keys():
        if 'Lum' in key:
            L += change[key]
    return L == 1



if __name__=="__main__":
     pass