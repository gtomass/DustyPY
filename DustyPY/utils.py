import matplotlib.pyplot as plt
import numpy as np
import subprocess
import glob
import pandas as pd
import os
import astropy.units as u
from astropy.table import Table, Column
from astropy.io import fits
from scipy.interpolate import make_interp_spline
import scipy.interpolate as interpolate
from pymcmcstat.MCMC import MCMC
from synphot import SpectralElement
import ctypes
import numpy as np
# ...existing imports...

# Load the compiled C library
simpson_lib = ctypes.CDLL('/Users/gtomassini/Documents/git_repo/DustyPY/DustyPY/libs/simpson.so')

# Define the argument and return types for the C functions
simpson_lib.simpson.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int]
simpson_lib.simpson.restype = ctypes.c_double

simpson_lib.simpson_error.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int]
simpson_lib.simpson_error.restype = ctypes.c_double

def simpson(x, y):
    """
    Python wrapper for the C implementation of Simpson's rule.

    Parameters:
    x (np.array): Array of x values.
    y (np.array): Array of y values.

    Returns:
    float: The integral calculated using Simpson's rule.
    """
    n = len(x) if len(x) % 2 == 1 else len(x) - 1
    x = np.array(x, dtype=np.double)
    y = np.array(y, dtype=np.double)

    x_c = x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    y_c = y.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    return simpson_lib.simpson(x_c, y_c, n)

def simpson_error(x, y):
    """
    Python wrapper for the C implementation of Simpson's rule error estimation.

    Parameters:
    x (np.array): Array of x values.
    y (np.array): Array of y values.

    Returns:
    float: The error estimate for the integral.
    """
    n = len(x) if len(x) % 2 == 1 else len(x) - 1
    x = np.array(x, dtype=np.double)
    y = np.array(y, dtype=np.double)

    x_c = x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    y_c = y.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    return simpson_lib.simpson_error(x_c, y_c, n)

def scatter_plot(Flux, Wavelength, unit=None, xlim=None, ylim=None, ax=None, scale='linear', kwargs=None, normalize: bool = False):
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
        unit = {'x': 'Wavelength', 'y': 'Flux'}
    if xlim is None:
        xlim = {}
    if ylim is None:
        ylim = {}
    if kwargs is None:
        kwargs = {}
    if ax is None:
        fig, ax = plt.subplots()

    if normalize:
        if xlim != {}:
            mask = (Wavelength >= xlim[0]) & (Wavelength <= xlim[1])
            Flux = Flux / np.max(Flux[mask])
        else:
            Flux = Flux / np.max(Flux)

    np.asarray(Wavelength) if type(Wavelength) == list else Wavelength

    ax.scatter(Wavelength, Flux, **kwargs)
    ax.set_xlabel(unit['x'], fontsize=15)
    ax.set_ylabel(unit['y'], fontsize=15)

    if scale == 'log-log':
        ax.set_xscale('log')
        ax.set_yscale('log')
    else:
        ax.set_yscale(scale)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)


def plot(Flux, Wavelength, unit=None, xlim=None, ylim=None, ax=None, scale='linear', kwargs=None, normalize: bool = False):
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
        unit = {'x': 'Wavelength', 'y': 'Flux'}
    if xlim is None:
        xlim = {}
    if ylim is None:
        ylim = {}
    if kwargs is None:
        kwargs = {}
    if ax is None:
        fig, ax = plt.subplots()

    if normalize:
        if xlim != {}:
            mask = (Wavelength >= xlim[0]) & (Wavelength <= xlim[1])
            Flux = Flux / np.max(Flux[mask])
        else:
            Flux = Flux / np.max(Flux)

    ax.plot(Wavelength, Flux, **kwargs)
    ax.set_xlabel(f"{unit['x']}",fontsize=15)
    ax.set_ylabel(f"{unit['y']}",fontsize=15)

    if scale == 'log-log':
        ax.set_xscale('log')
        ax.set_yscale('log')
    else:
        ax.set_yscale(scale)


    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)


def error_plot(Flux, Wavelength, eFlux, unit=None, xlim=None, ylim=None, ax=None, scale='linear', kwargs=None, normalize: bool = False):
    """Create a scatter plot
    """

    if unit is None:
        unit = {'x': 'Wavelength', 'y': 'Flux'}
    if xlim is None:
        xlim = {}
    if ylim is None:
        ylim = {}
    if kwargs is None:
        kwargs = {}
    if ax is None:
        fig, ax = plt.subplots()

    if normalize:
        if xlim != {}:
            mask = (Wavelength >= xlim[0]) & (Wavelength <= xlim[1])
            Flux = Flux / np.max(Flux[mask])
        else:
            Flux = Flux / np.max(Flux)

    ax.errorbar(Wavelength, Flux, yerr=eFlux, **kwargs)
    ax.set_xlabel(f"{unit['x']}",fontsize=15)
    ax.set_ylabel(f"{unit['y']}",fontsize=15)

    if scale == 'log-log':
        ax.set_xscale('log')
        ax.set_yscale('log')
    else:
        ax.set_yscale(scale)


    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)


def print_file(file, stop=-1):
    """Print a file

    Args:
        file (string): Path to the file
        stop (int, optional): line where to stop printing. Defaults to -1.
    """
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines[:stop]:
            print(line)


def load_file(Path, header: int = 0):
    """Load a file in an array

    Args:
        Path (string): Path to the file 

    Returns:
        array: array containing the lines of the file
    """

    with open(Path, 'r') as f:
        return f.readlines()[header:]


def load_csv(Path, sep=','):
    """Load a csv file

    Args:
        Path (string): Path to the file 

    Returns:
        Dataframe: dataframe containing the csv file
    """
    return pd.read_csv(Path, sep=sep)


def load_fits(Path):
    """Load a fits file

    Args:
        Path (string): Path to the file 

    Returns:
        array: array containing the file
    """
    return fits.open(Path)[0].data


def search_line(file, line):
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
    for i, lines in enumerate(file):
        if line in lines:
            return i
    raise Exception(f'This line does dot exist: {line}')


def save_file(Path, file):
    """
    Save the given content to a file at the specified path.

    Args:
        Path (str): The path where the file will be saved.
        file (list of str): The content to be written to the file, provided as a list of strings.

    Returns:
        None
    """
    with open(Path, 'w') as f:
        f.write("".join(file))


def build_change_dict(model):
    """
    Build a dictionary with various parameters extracted from the model.

    Args:
        model: The model object containing the necessary data.

    Returns:
        dict: A dictionary with the extracted parameters.
    """
    stars = model.get_Stars()
    dust = model.get_Dust()
    composition = dust.get_Composition()

    temperatures = ", ".join(str(star.get_Temperature()) for star in stars)
    luminosities = ", ".join(str(star.get_Luminosity()) for star in stars)
    composition_files = "\n        ".join(
        f"{os.path.join('data', 'Lib_nk', comp)}.nk" for comp in composition.keys())
    abundances = ", ".join(str(ab) for ab in composition.values())
    
    if dust.get_Density()['density type'] == 'POWD':
        density = f"\t\t density type = {dust.get_Density()['density type']}\n \t\t number of powers = {dust.get_Density()['number of powers']}\n \t\t shell's relative thickness = {dust.get_Density()['shell']}\n \t\t power = {dust.get_Density()['power']}\n"
    elif dust.get_Density()['density type'] == 'RDWA':
        density = f"\t\t density type = {dust.get_Density()['density type']} ;\n \t\t Y = {dust.get_Density()['shell']}\n"
    elif dust.get_Density()['density type'] == 'RDW':
        density = f"\t\t density type = {dust.get_Density()['density type']} ;\n \t\t  Y = {dust.get_Density()['shell']}\n"
    else:
        raise NotImplementedError(f'This density type is not implemented: {dust.get_Density()["density type"]}')

    return {
        'Spectral': f'      	        Spectral shape = {model.get_Spectral()} \n' if model.get_Spectral() in ['black_body', 'engelke_marengo'] else
        f'      	        Spectral shape = {model.get_Spectral()} \n \t\t{
            model.get_SpectralFile()} \n',
        'BB': f'        	Number of BB = {len(stars)} \n',
        'Temperature': f'        	Temperature = {temperatures} K \n',
        'Luminosities': f'        	Luminosities = {luminosities} \n',
        'Dust Temperature': f'        		Scale:    type of entry = T1\n  \t\t\t Td = {dust.get_Temperature()} K \n',
        'Absorption': f'        SiO absorption depth = {model.get_SiOAbsorption()}  percents\n',
        'Optical properties': f'        optical properties index = {dust.get_Properties()} \n',
        'Composition': f'	Number of additional components = {len(composition)}, properties listed files \n        {composition_files}\n',
        'Abundances': f'   Abundances for these components = {abundances} \n',
        'Size Distribution': f'        SIZE DISTRIBUTION = {dust.get_DustSize()["Distribution"]} \n',
        'Dust size': f'        q = {dust.get_DustSize()["q"]}, a(min) = {dust.get_DustSize()["amin"]} micron, a(max) = {dust.get_DustSize()["amax"]} micron \n',
        'Sublimation temperature': f'        Tsub = {dust.get_Sublimation()} K \n',
        'Density Distribution': density,
        'Opacity': f'        - tau(min) = {dust.get_tau()}; tau(max) = {dust.get_tau()}  % for the visual wavelength \n',
    }


def change_parameter(Path, change, car, nstar):
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
    file = load_file(Path)

    # Remove all lines containing '.nk'
    file = [line for line in file if ('.nk' not in line)]
    file = [line for line in file if ('Td' not in line)]
    file = [line for line in file if all(param not in line for param in ['number of powers', 'power =', 'shell', 'Y ='])]

    # cannot have multiple stars with engelke_marengo
    if 'engelke_marengo' in change['Spectral']:
        change.pop('BB')
        change.pop('Luminosities')
        file = [line for line in file if ('Number of BB' not in line)]

    elif 'black_body' in change['Spectral']:
        change.pop('Absorption')
        file = [line for line in file if ('SiO absorption depth' not in line)]

    else:
        change.pop('BB')
        change.pop('Absorption')
        change.pop('Luminosities')
        change.pop('Temperature')
        file = [line for line in file if all(param not in line for param in [
                                             'Number of BB', 'SiO absorption depth', 'Luminosities', 'Temperature'])]

    if ('MRN' in change['Size Distribution']) and ('MODIFIED_MRN' not in change['Size Distribution']):
        change.pop('Dust size')
        file = [line for line in file if ('q = 3.5' not in line)]

    if nstar == 1 and 'Luminosities' in change.keys():
        file = [line for line in file if ('Luminosities' not in line)]
        change.pop('Luminosities')

    for param in change.keys():
        line = search_line(file, car[param])
        new_line = change[param]
        file[line] = new_line

    save_file(Path, file)

def set_change(dusty, change: dict) -> None:
    """
    Applies changes to the Dusty model based on the provided dictionary.

    Parameters:
    change (dict): A dictionary containing the changes to be applied. The keys should specify the parameter to change (e.g., 'Temp', 'Lum', 'Opacity' for each stars of the model) and the values should be the new values for those parameters.

    Raises:
    NotImplementedError: If an attempt is made to change the dust size, which is not yet fittable.
    """
    for key in change.keys():
        if 'Temp' in key and key != 'Temperature':
            dusty.get_Model().get_Stars()[int(
                key.split('Temp')[-1])-1].set_Temperature(change[key])
        elif 'Lum' in key:
            dusty.get_Model().get_Stars()[int(
                key.split('Lum')[-1])-1].set_Luminosity(change[key])
        elif key in ['Opacity']:
            dusty.get_Model().get_Dust().set_tau(change[key])
        elif key in ['Composition', 'Abundances']:
            # self._Dusty.get_Model().get_Dust().set_Composition(change[key]
            # Composition must be fixed but abundances can be fitted
            raise NotImplementedError(
                'Composition and Abundances are not yet fittable.')
        elif key in ['DustSize']:
            dusty.get_Model().get_Dust().set_DustSize(change[key])
        elif key in ['Sublimation']:
            dusty.get_Model().get_Dust().set_Sublimation(
                change[key])
        elif key in ['Absorption']:
            dusty.get_Model().set_SiOAbsorption(change[key])
        elif key in ['Temperature']:
            dusty.get_Model().get_Dust().set_Temperature(change[key])
        elif key in ['Density']:
            dusty.get_Model().get_Dust().set_Density(change[key])

        elif key == 'Lest':
            pass
        else:
            raise NotImplementedError(f'Parameter {key} not recognized.')


def list_to_dict(keys, values):
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
    >>> list_to_dict(keys, values)
    {'a': 1, 'b': 2, 'c': 3}
    """
    key_value_pairs = zip(keys, values)
    return dict(key_value_pairs)


def str_to_data(data=list) -> np.array:
    """
    Convert a list of strings into a numpy array.

    This function takes a list of strings and converts it into a numpy array.

    Parameters:
    data (list): A list of strings.

    Returns:
    np.array: A numpy array containing the data from the input list.

    Example:
    >>> data = ['1 2 3', '4 5 6', '7 8 9']
    >>> str_to_data(data)
    array([[1, 2, 3],
           [4, 5, 6],
           [7, 8, 9]])
    """
    return np.array([list(map(float, line.split())) for line in data])


def get_column_spectum(file, index, index_header=0):
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
    return np.array([[x for x in el.split('  ')[1:-1] if x != ''] for el in array], dtype=float).T[index]


def watt_to_jansky(Flux, Wavelength):
    """
    Convertit le flux de Watts par mètre carré en Jansky.

    Paramètres:
    Flux (float): Le flux en Watts par mètre carré.
    Wavelength (float): La longueur d'onde en micromètres.

    Retourne:
    float: Le flux en Jansky.
    """
    return (Flux * u.W/u.m**2).to(u.Jy, equivalencies=u.spectral_density(Wavelength * u.um)).value


def jansky_to_watt(Flux, Wavelength):
    """
    Convertit le flux de Jansky en Watts par mètre carré.

    Paramètres:
    Flux (float): Le flux en Jansky.
    Wavelength (float): La longueur d'onde en micromètres.

    Retourne:
    float: Le flux en Watts par mètre carré.
    """
    return (Flux * u.Jy).to(u.W/u.m**2, equivalencies=u.spectral_density(Wavelength * u.um)).value


def um_to_meter(Wavelength):
    """
    Convertit la longueur d'onde de micromètres en mètres.

    Paramètres:
    Wavelength (float): La longueur d'onde en micromètres.

    Retourne:
    float: La longueur d'onde en mètres.
    """
    return (Wavelength * u.um).to(u.m).value


def calcul_rayon_vrai(results, L):
    """
    Calcule le rayon réel en fonction des résultats et de la luminosité.

    Paramètres:
    results (dict): Un dictionnaire contenant les résultats, y compris 'r1(cm)'.
    L (float): La luminosité en unités arbitraires.

    Retourne:
    float: Le rayon réel en mètres.
    """
    return results['r1(cm)'] * np.sqrt(L / 1e4) * 1e-2


def calcul_flux_total(F, r_vrai, distance):
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


def interpolate_spline(Wavelength, Flux, order=3):
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


def interpol(xdata, xdusty, ydusty):
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
    try:
        return interpolate_spline(np.asarray(xdusty).flatten(), np.asarray(ydusty).flatten())(np.asarray(xdata).flatten())
    except Exception:
        return interpolate_spline(xdusty, ydusty)(xdusty)
    
def model(theta, data)-> None:
        
        dusty, data_mod, fit, logfile, Jansky = data.user_defined_object[0]

        dustsize = dusty.get_Model().get_Dust().get_DustSize()
        density = dusty.get_Model().get_Dust().get_Density()
        p = dict(zip(fit.get_Param().keys(), theta))

        if 'amin' in p:
            dustsize['amin'] = np.round(10**p['amin'], 3)
        if 'amax' in p:
            dustsize['amax'] =  np.round(10**p['amax'], 3)
        if 'q' in p:
            dustsize['q'] = p['q']
        if 'shell' in p:
            density['shell'] = p['shell']


        if dustsize['amin'] > dustsize['amax']:
            ymodel = np.zeros_like(data.xdata)
            ymodel[:] = np.nan
            return ymodel
        
        if density['shell'] < 0:
            ymodel = np.zeros_like(data.xdata)
            ymodel[:] = np.nan
            return ymodel

        change = {key: value for key, value in list_to_dict(list(fit.get_Param().keys()), theta).items() if key not in ['amin', 'amax', 'q', 'shell']}

        change.update({'DustSize': dustsize})
        change.update({'Density': density})
        Lum  = p['Lest']

        set_change(dusty,change)

        dusty.change_parameter()
        dusty.lunch_dusty(verbose=0, logfile=logfile)
        dusty.make_SED(luminosity=Lum, Jansky=Jansky)

        if data_mod.get_table() is not None:

            bandpass = data_mod.get_common_filters(data_mod.get_table())
            
            ymodel,ymodel_err  =  dusty.get_SED().integrate_bandpass(bandpass=bandpass)

            central_wavelength = [get_central_wavelegnth(get_bandpass(f)) for f in bandpass.values()]
            index = np.argsort(central_wavelength)
            ymodel = ymodel[index]

        else:
            ymodel = interpol(data.xdata, dusty.get_SED().get_Wavelength(), dusty.get_SED().get_Flux()
                             )
            
        subprocess.call('clear', shell=True)

        return ymodel

def prediction_model(theta, data):
    """
    Calcule le modèle en fonction des paramètres fournis.

    Paramètres:
    theta (array-like): Les paramètres du modèle.
    data (tuple or object): Les données observées. Peut être un tuple (xdata, ydata) ou un objet avec des attributs xdata et ydata.

    Retourne:
    array-like: Les valeurs du modèle interpolé.
    """

    dusty, data_mod, fit, logfile, Jansky = data.user_defined_object[0]

    dustsize = dusty.get_Model().get_Dust().get_DustSize()
    density = dusty.get_Model().get_Dust().get_Density()
    p = dict(zip(fit.get_Param().keys(), theta))

    if 'amin' in p:
        dustsize['amin'] = np.round(10**p['amin'], 3)
    if 'amax' in p:
        dustsize['amax'] =  np.round(10**p['amax'], 3)
    if 'q' in p:
        dustsize['q'] = p['q']
    if 'shell' in p:
        density['shell'] = p['shell']

    if dustsize['amin'] > dustsize['amax']:
        ymodel = np.zeros_like(data.xdata)
        ymodel[:] = np.nan
        return ymodel
    
    if density['shell'] < 0:
        ymodel = np.zeros_like(data.xdata)
        ymodel[:] = np.nan
        return ymodel


    change = {key: value for key, value in list_to_dict(list(fit.get_Param().keys()), theta).items() if key not in ['amin', 'amax', 'q', 'shell']}
    change.update({'DustSize': dustsize})
    change.update({'Density': density})
    Lum  =  p['Lest']

    set_change(dusty,change)

    dusty.change_parameter()
    dusty.lunch_dusty(verbose=0, logfile=logfile)
    dusty.make_SED(luminosity=Lum, Jansky=Jansky)
    subprocess.call('clear', shell=True)
    ymodel = interpol(data.xdata, dusty.get_SED().get_Wavelength(), dusty.get_SED().get_Flux())
    return ymodel


def chi2(theta, data):
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

    ymodel = interpol(theta, xdata, ydata).reshape(ydata.shape)
    return np.nansum((ymodel - ydata)**2)


def set_mcmc_param(mc=MCMC, param=None):
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
                                          maximum=param[par]['maximum'],
                                          sample=param[par]['sample'])
        
def supp_car_list(list, car):
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


def unred(Wavelength, Flux, EBV, Rv=3.1, LMC2=False, AVGLMC=False) -> np.array:
    """
    Corrects flux for interstellar extinction using Fitzpatrick (1999) extinction curve.

    Parameters:
    Wavelength : array-like
        Wavelengths in Angstroms.
    Flux : array-like
        Flux values corresponding to the input wavelengths.
    EBV : float
        Color excess E(B-V).
    Rv : float, optional
        Total-to-selective extinction ratio (default 3.1).
    LMC2 : bool, optional
        Use LMC2 extinction curve (default False).
    AVGLMC : bool, optional
        Use average LMC extinction curve (default False).

    Returns:
    array-like
        Flux values corrected for interstellar extinction.
    """
    x = 10000. / Wavelength  # Convert to inverse microns
    curve = np.zeros_like(x)

    # Default parameters
    params = {
        "x0": 4.596, "gamma": 0.99, "c3": 3.23, "c4": 0.41,
        "c2": -0.824 + 4.717 / Rv, "c1": 2.030 - 3.007 * (-0.824 + 4.717 / Rv)
    }

    if LMC2:
        params.update({"x0": 4.626, "gamma": 1.05, "c3": 1.92, "c4": 0.42, "c2": 1.31, "c1": -2.16})
    elif AVGLMC:
        params.update({"x0": 4.596, "gamma": 0.91, "c3": 2.73, "c4": 0.64, "c2": 1.11, "c1": -1.28})

    # UV portion
    xcutuv = 10000.0 / 2700.0
    iuv = np.where(x >= xcutuv)[0]
    xuv = np.concatenate(([10000.0 / 2700.0, 10000.0 / 2600.0], x[iuv])) if len(iuv) > 0 else [10000.0 / 2700.0, 10000.0 / 2600.0]
    xuv = np.array(xuv)  # Ensure xuv is a NumPy array
    yuv = params["c1"] + params["c2"] * xuv + params["c3"] * xuv**2 / ((xuv**2 - params["x0"]**2)**2 + (xuv * params["gamma"])**2)
    yuv += params["c4"] * (0.5392 * (np.maximum(xuv, 5.9) - 5.9)**2 + 0.05644 * (np.maximum(xuv, 5.9) - 5.9)**3) + Rv
    if len(iuv) > 0:
        curve[iuv] = yuv[2:]

    # Optical/IR portion
    xsplopir = np.concatenate(([0], 10000.0 / np.array([26500.0, 12200.0, 6000.0, 5470.0, 4670.0, 4110.0])))
    ysplopir = np.concatenate((np.array([0.0, 0.26469, 0.82925]) * Rv / 3.1, [np.polyval([-0.422809, 1.0027, 0.000213572][::-1], Rv),
                                                                             np.polyval([-0.051354, 1.00216, -0.000073578][::-1], Rv),
                                                                             np.polyval([0.700127, 1.00184, -0.00003326][::-1], Rv),
                                                                             np.polyval([1.19456, 1.01707, -0.00546959, 0.000797809, -0.0000445636][::-1], Rv)]))
    iopir = np.where(x < xcutuv)[0]
    if len(iopir) > 0:
        tck = interpolate.splrep(np.concatenate((xsplopir, [10000.0 / 2700.0, 10000.0 / 2600.0])), np.concatenate((ysplopir, yuv[:2])), s=0)
        curve[iopir] = interpolate.splev(x[iopir], tck)

    return Flux * 10.**(0.4 * curve * EBV)


def querry_vizier_data(radius, target):
    """
    Interroge Vizier pour obtenir des données.

    Paramètres:
    radius (str): Le rayon de la requête.
    target (str): La cible de la requête.

    Retourne:
    array-like: Les données obtenues.
    """

    return Table.read(f"https://vizier.cds.unistra.fr/viz-bin/sed?-c={target}&-c.rs={radius}")


def log_space(start, stop, num):
    """
    Crée un tableau de valeurs espacées logarithmiquement.

    Paramètres:
    start (float): La valeur de départ.
    stop (float): La valeur de fin.
    num (int): Le nombre de valeurs.

    Retourne:
    array-like: Un tableau de valeurs espacées logarithmiquement.
    """
    return np.logspace(np.log10(start), np.log10(stop), num)


def write_wavelength(Path, Wavelength):
    """
    Écrit les longueurs d'onde dans un fichier.

    Paramètres:
    Path (str): Le chemin du fichier.
    Wavelength (array-like): Les longueurs d'onde à écrire.

    Retourne:
    None
    """
    header = []

    with open(Path, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if '# nL =' in line:
                line = f'# nL = {len(Wavelength)}\n'
                header.append(line)
            elif '#' in line:
                header.append(line)
                pass
            else:
                break

    with open(Path, 'w') as f:
        for line in header:
            f.write(line)
        for w in Wavelength:
            f.write(f"{w}\n")

def get_bandpass_name() -> list:
    """
    Get the name of all the available filters.

    Returns:
    list: A list of filter names.
    """
    path_here = os.path.dirname(__file__)
    test = glob.glob(os.path.join(path_here, 'filter', 'comp', 'nonhst', '*'))
    return [os.path.basename(f).split('.')[0] for f in test]

def get_bandpass(bandpass: str) -> SpectralElement:
    """
    Get the bandpass of a filter.

    Parameters:
    filter (str): The name of the filter to retrieve the bandpass for.

    Returns:
    SpectralElement: The spectral element corresponding to the filter.

    Raises:
    FileNotFoundError: If the specified filter is not found in the directory.
    """
    path_here = os.path.dirname(__file__)
    test = glob.glob(os.path.join(path_here, 'filter', 'comp', 'nonhst', '*'))
    filter_path = next((f for f in test if f.endswith(f'{bandpass}.fits')), None)

    if filter_path is None:
        raise FileNotFoundError(f"Filter {bandpass} not found")
    return SpectralElement.from_file(filter_path)

def intergrate_bandpass(wavelength: np.array, flux: np.array, bandpass: SpectralElement, n_int: int = 10000) -> float:
    """
    Integrate the bandpass of a given spectrum using the custom Simpson's rule.

    Parameters:
    wavelength (np.array): Array of wavelength values.
    flux (np.array): Array of flux values corresponding to the wavelengths.
    bandpass (SpectralElement): Filter object that defines the bandpass.

    Returns:
    float: The integrated flux over the bandpass.
    """
    wavelength_int = np.linspace(np.min(bandpass.waveset.value), np.max(bandpass.waveset.value), n_int) *u.AA
    flux_interp = np.interp(wavelength_int.value / 10000, wavelength, flux)
    filtSpec = bandpass(wavelength_int) * flux_interp  # Calculate throughput
    numerator, numerator_err = simpson(wavelength_int.value / 10000, filtSpec), simpson_error(wavelength_int.value / 10000, bandpass(wavelength_int))
    denominator, denominator_err = simpson(wavelength_int.value / 10000, bandpass(wavelength_int)), simpson_error(wavelength_int.value / 10000, bandpass(wavelength_int))
    if denominator == 0 :
        fl,fl_err2 = 0, 0
    elif denominator_err == 0:
        fl = numerator / denominator
        fl_err2 = 0
    else:
        fl = numerator / denominator
        fl_err2 = (1/denominator)**2*numerator_err**2 + (numerator/denominator**2)**2*denominator_err**2
    return fl, np.sqrt(fl_err2)

def get_central_wavelegnth(bandpass: SpectralElement) -> float:
    """
    Get the central wavelength of a bandpass.

    Parameters:
    bandpass (SpectralElement): The bandpass object.

    Returns:
    float: The central wavelength of the bandpass.
    """
    return bandpass.pivot().value 

def mean_flux(wavelength: np.array, flux: np.array) -> float:
    """
    Calculate the mean flux value over a specified wavelength range.

    Parameters:
    wavelength (np.array): Array of wavelength values.
    flux (np.array): Array of flux values corresponding to the wavelengths.

    Returns:
    float: The mean flux value over the specified wavelength range.
    """
    unique_wavelengths = np.unique(wavelength)
    mean_flux_values = []

    for wl in unique_wavelengths:
        mask = (wavelength == wl)
        mean_flux_values.append(np.mean(flux[mask]))

    return np.asarray(mean_flux_values)

def get_common_filters(filters: np.array, bandpass_name: np.array) -> dict:
    """
    Get the common filters between the filters and bandpass names.

    Parameters:
    filters (np.array): Array of filter names to check.
    bandpass_name (np.array): Array of available bandpass names.

    Returns:
    dict: A dictionary with filter names as keys and corresponding bandpass names as values.
    """
    bandpass_name = ['_'.join(name.split('_')[:2]) for name in bandpass_name]
    common_filters = {}
    for f in filters:
        for bp in bandpass_name:
            if bp in f.lower().replace(':', '_') and bp not in common_filters.values():
                common_filters[f] = bp
                break
    return common_filters

def integrate_SED_bandpass(wavelength: np.array, flux: np.array, common_filter: np.array) -> np.array:
    """
    Integrate the bandpass of a given spectrum.

    Parameters:
    wavelength (np.array): Array of wavelength values.
    flux (np.array): Array of flux values corresponding to the wavelengths.
    filter (np.array): Array of filter values corresponding to the wavelengths.

    Returns:
    np.array: The integrated flux over the bandpass.
    """
    integrated_flux, integrated_flux_err = [], []

    for bandpass in common_filter.values():
        filt  = get_bandpass(bandpass)
        fl, fl_err = intergrate_bandpass(wavelength, flux, filt)
        integrated_flux.append(fl)
        integrated_flux_err.append(fl_err)
        
    
    return np.asarray(integrated_flux) , np.asarray(integrated_flux_err) 

def savefig(path: str) -> None: 
    """
    Save the current figure to a file.

    Parameters:
    path (str): The path to save the figure to.

    Returns:
    None
    """
    plt.savefig(path) 

def write_table_to_latex(table, Path, columns=None, column_names=None, wavelength: bool = True) -> None:
    """
    Write a table to a LaTeX file.

    Parameters:
    table (DataFrame): The table to write to the file.
    Path (str): The path to save the table to.
    columns (list, optional): List of columns to include in the LaTeX file. Defaults to None (all columns).
    column_names (list, optional): List of column names to use in the LaTeX file. Defaults to None (original names).

    Returns:
    None
    """
    if wavelength:
        table.add_column(Column(data=table['sed_freq'].to(u.um, equivalencies=u.spectral()), name='Wavelength'),)
    if columns is None:
        columns = table.columns
    if column_names is None:
        column_names = columns

    with open(Path, 'w') as f:
        f.write('\\begin{table*}[h!]\n')
        f.write('\\centering\n')
        f.write('\\begin{tabular*}{\\textwidth}{@{\extracolsep{\\fill}}' + 'c' * len(columns) + '}\n')
        f.write('\\hline\n')
        f.write(' & '.join(column_names) + ' \\\\\n')
        f.write('\\hline\n')
        for row in table[columns]:
            f.write(' & '.join(f'{value:.3f}' if isinstance(value, (int, float)) and value != 0 else '-' if value == 0 else str(value) for value in row) + ' \\\\\n')
        f.write('\\hline\n')
        f.write('\\end{tabular*}\n')
        f.write('\\caption{Your caption here}\n')
        f.write('\\label{table:label}\n')
        f.write('\\end{table*}\n')

def get_photometry(Wavelength, Flux, bandpass_name: str) -> tuple:
    """
    Get the photometry of a given spectrum.

    Parameters:
    Wavelength (np.array): Array of wavelength values.
    Flux (np.array): Array of flux values corresponding to the wavelengths.
    bandpass_name (str): The name of the bandpass to use.

    Returns:
    tuple: A tuple containing the integrated flux and its error.
    """
    bandpass = get_bandpass(bandpass_name)
    integrated_flux, integrated_flux_err = intergrate_bandpass(Wavelength, Flux, bandpass)    
    return get_central_wavelegnth(bandpass)/10000,integrated_flux, integrated_flux_err



if __name__ == "__main__":
    pass
