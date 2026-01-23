import os
import glob
import ctypes
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.table import Table, Column
from astropy.io import fits
from scipy.interpolate import make_interp_spline, LinearNDInterpolator
import scipy.interpolate as interpolate
from synphot import SpectralElement
from typing import Callable, Dict, List, Tuple, Any, Optional

# Import Prior from your new PyMCMC library
from PyMCMC import Prior

# =============================================================================
# C LIBRARY LOADING (Simpson's Rule)
# =============================================================================

def load_simpson_library():
    """
    Dynamically loads the compiled simpson shared library (.so file).
    """
    lib_dir = os.path.join(os.path.dirname(__file__), "libs")
    for file in os.listdir(lib_dir):
        if file.startswith("simpson") and file.endswith(".so"):
            return ctypes.CDLL(os.path.join(lib_dir, file))
    raise ImportError("The simpson shared library could not be found. Ensure it is compiled correctly.")

try:
    simpson_lib = load_simpson_library()
    simpson_lib.simpson.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int]
    simpson_lib.simpson.restype = ctypes.c_double
    simpson_lib.simpson_error.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int]
    simpson_lib.simpson_error.restype = ctypes.c_double
except ImportError as e:
    raise ImportError("The simpson shared library could not be loaded.") from e

def simpson(x: np.ndarray, y: np.ndarray) -> float:
    """Python wrapper for the C implementation of Simpson's rule."""
    n = len(x) if len(x) % 2 == 1 else len(x) - 1
    x_c = np.array(x, dtype=np.double).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    y_c = np.array(y, dtype=np.double).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    return simpson_lib.simpson(x_c, y_c, n)

def simpson_error(x: np.ndarray, y: np.ndarray) -> float:
    """Python wrapper for the C implementation of Simpson's rule error estimation."""
    n = len(x) if len(x) % 2 == 1 else len(x) - 1
    x_c = np.array(x, dtype=np.double).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    y_c = np.array(y, dtype=np.double).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    return simpson_lib.simpson_error(x_c, y_c, n)

# =============================================================================
# PARAMETER CONVERSION
# =============================================================================

def convert_params_to_priors(param_dict: Dict[str, Dict[str, Any]]) -> Tuple[List[Prior], List[float]]:
    """
    Converts Dusty parameter dictionary to PyMCMC Prior objects.

    Args:
        param_dict: Dictionary defining parameters with 'sample', 'minimum', 
                    'maximum', and 'theta0' keys.

    Returns:
        A list of Prior objects and a list of initial values for sampled parameters.
    """
    priors = []
    initial_theta0 = []
    for name, settings in param_dict.items():
        if settings.get('sample', False):
            # Standard Dusty limits are treated as Uniform priors
            if settings.get('prior_type', 'uniform') == 'uniform':
                priors.append(Prior('uniform', (settings['minimum'], settings['maximum'])))
                initial_theta0.append(settings['theta0'])
            elif settings['prior_type'] == 'gaussian':
                priors.append(Prior('gaussian', (settings['mu'], settings['sigma'])))
                initial_theta0.append(settings['theta0'])
            elif settings['prior_type'] == 'log_normal':
                priors.append(Prior('log_normal', (settings['mu'], settings['sigma'])))
                initial_theta0.append(settings['theta0'])
            elif settings['prior_type'] == 'beta':
                priors.append(Prior('beta', (settings['a'], settings['b'])))
                initial_theta0.append(settings['theta0'])
            else:
                raise ValueError(f"Unsupported prior type: {settings['prior_type']}")

    return priors, initial_theta0

# =============================================================================
# PLOTTING UTILITIES
# =============================================================================

def scatter_plot(Flux, Wavelength, unit=None, xlim=None, ylim=None, ax=None, scale='linear', kwargs=None, normalize: bool = False):
    """Creates a scatter plot for flux vs wavelength."""
    if unit is None: unit = {'x': 'Wavelength ($\mu m$)', 'y': 'Flux'}
    if kwargs is None: kwargs = {}
    if ax is None: fig, ax = plt.subplots()

    if normalize:
        if xlim:
            mask = (Wavelength >= xlim[0]) & (Wavelength <= xlim[1])
            Flux = Flux / np.max(Flux[mask])
        else:
            Flux = Flux / np.max(Flux)

    ax.scatter(Wavelength, Flux, **kwargs)
    ax.set_xlabel(unit['x'], fontsize=15)
    ax.set_ylabel(unit['y'], fontsize=15)
    if scale == 'log-log':
        ax.set_xscale('log')
        ax.set_yscale('log')
    else:
        ax.set_yscale(scale)
    if xlim: ax.set_xlim(*xlim)
    if ylim: ax.set_ylim(*ylim)

def plot(Flux, Wavelength, unit=None, xlim=None, ylim=None, ax=None, scale='linear', kwargs=None, normalize: bool = False):
    """Creates a standard line plot for flux vs wavelength."""
    if unit is None: unit = {'x': 'Wavelength ($\mu m$)', 'y': 'Flux'}
    if kwargs is None: kwargs = {}
    if ax is None: fig, ax = plt.subplots()

    if normalize:
        if xlim:
            mask = (Wavelength >= xlim[0]) & (Wavelength <= xlim[1])
            Flux = Flux / np.max(Flux[mask])
        else:
            Flux = Flux / np.max(Flux)

    ax.plot(Wavelength, Flux, **kwargs)
    ax.set_xlabel(unit['x'], fontsize=15)
    ax.set_ylabel(unit['y'], fontsize=15)
    if scale == 'log-log':
        ax.set_xscale('log')
        ax.set_yscale('log')
    else:
        ax.set_yscale(scale)
    if xlim: ax.set_xlim(*xlim)
    if ylim: ax.set_ylim(*ylim)

def error_plot(Flux, Wavelength, eFlux, unit=None, xlim=None, ylim=None, ax=None, scale='linear', kwargs=None, normalize: bool = False):
    """Creates a plot with error bars for flux vs wavelength."""
    if unit is None: unit = {'x': 'Wavelength ($\mu m$)', 'y': 'Flux'}
    if kwargs is None: kwargs = {}
    if ax is None: fig, ax = plt.subplots()

    if normalize:
        if xlim:
            mask = (Wavelength >= xlim[0]) & (Wavelength <= xlim[1])
            Flux = Flux / np.max(Flux[mask])
        else:
            Flux = Flux / np.max(Flux)

    ax.errorbar(Wavelength, Flux, yerr=eFlux, **kwargs)
    ax.set_xlabel(unit['x'], fontsize=15)
    ax.set_ylabel(unit['y'], fontsize=15)
    if scale == 'log-log':
        ax.set_xscale('log')
        ax.set_yscale('log')
    else:
        ax.set_yscale(scale)
    if xlim: ax.set_xlim(*xlim)
    if ylim: ax.set_ylim(*ylim)

# =============================================================================
# FILE HANDLING
# =============================================================================

def load_file(Path: str, header: int = 0) -> List[str]:
    """Loads a text file into a list of strings."""
    with open(Path, 'r') as f:
        return f.readlines()[header:]

def load_csv(Path: str, sep: str = ',') -> pd.DataFrame:
    """Loads a CSV file into a pandas DataFrame."""
    return pd.read_csv(Path, sep=sep)

def load_fits(Path: str) -> np.ndarray:
    """Loads a FITS file and returns the data from the primary HDU."""
    return fits.open(Path)[0].data

def search_line(file: List[str], line_content: str) -> int:
    """Searches for a specific substring in a list of lines and returns the index."""
    for i, line in enumerate(file):
        if line_content in line:
            return i
    raise ValueError(f"Line containing '{line_content}' not found.")

def get_table_interpolated(teff=None, logg=None, ebv=0.0, **kwargs) -> tuple:
    """
    Interpolates the atmosphere model grid to the desired Teff and logg.
    """
    # -- get the grid

    gridfilename = kwargs.get('gridfilename', None)
    hdu = fits.open(gridfilename)

    wave_unit = hdu[0].header['WAVUNIT']
    if wave_unit.lower() == 'angstrom':
        wave_unit = u.Angstrom
    elif wave_unit.lower() == 'micron':
        wave_unit = u.micron
    elif wave_unit.lower() == 'nm':
        wave_unit = u.nm
    else:
        raise ValueError(f"Unknown wavelength unit: {wave_unit}")

    # Extract grid points and models
    teffs, loggs, models = [], [], []
    for i in range(1, len(hdu)):
        teffs.append(hdu[i].header['TEFF'])
        loggs.append(hdu[i].header['LOGG'])
        models.append(hdu[i].data)

    teffs = np.array(teffs)
    loggs = np.array(loggs)

    # Interpolate flux for each wavelength
    wave = (models[0]['wavelength'] * wave_unit).to(u.micron).value
    flux_grid = np.array([model['flux'] for model in models])
    # Find the closest values for teff and logg in the grid
    interpolator = LinearNDInterpolator(list(zip(teffs, loggs)), flux_grid)

    # Interpolate to the desired Teff and logg
    flux = interpolator(teff, logg)

    if flux is None:
        return False, wave, flux 

    return True, wave, flux


def save_file(Path: str, file_lines: List[str]) -> None:
    """Saves a list of strings to a file."""
    with open(Path, 'w') as f:
        f.write("".join(file_lines))

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

# =============================================================================
# DUSTY PARAMETER LOGIC
# =============================================================================

def build_change_dict(model) -> Dict[str, str]:
    """Builds a dictionary of formatted strings representing parameters for Dusty input files."""
    stars = model.get_Stars()
    dust = model.get_Dust()
    composition = dust.get_Composition()
    temperatures = ", ".join(str(star.get_Temperature()) for star in stars)
    luminosities = ", ".join(str(star.get_Luminosity()) for star in stars)
    composition_files = "\n        ".join(f"{os.path.join('data', 'Lib_nk', comp)}.nk" for comp in composition.keys())
    abundances = ", ".join(str(ab) for ab in composition.values())
    
    density_settings = dust.get_Density()
    if density_settings['density type'] == 'POWD':
        density = f"\t\t density type = POWD\n \t\t number of powers = {density_settings['number of powers']}\n \t\t shell's relative thickness = {density_settings['shell']}\n \t\t power = {density_settings['power']}\n"
    elif density_settings['density type'] in ['RDWA', 'RDW']:
        density = f"\t\t density type = {density_settings['density type']} ;\n \t\t Y = {density_settings['shell']}\n"
    else:
        raise NotImplementedError(f"Density type '{density_settings['density type']}' not implemented.")

    return {
        'Spectral': f"      Spectral shape = {model.get_Spectral()} \n" if model.get_Spectral() in ['black_body', 'engelke_marengo'] else f"      Spectral shape = {model.get_Spectral()} \n \t\t{model.get_SpectralFile()} \n",
        'BB': f"        Number of BB = {len(stars)} \n",
        'Temperature': f"        Temperature = {temperatures} K \n",
        'Luminosities': f"        Luminosities = {luminosities} \n",
        'Dust Temperature': f"        Scale: type of entry = T1\n \t\t Td = {dust.get_Temperature()} K \n",
        'Absorption': f"        SiO absorption depth = {model.get_SiOAbsorption()} percents\n",
        'Optical properties': f"        optical properties index = {dust.get_Properties()} \n",
        'Composition': f"        Number of additional components = {len(composition)}, properties listed files \n {composition_files}\n",
        'Abundances': f"        Abundances for these components = {abundances} \n",
        'Size Distribution': f"        SIZE DISTRIBUTION = {dust.get_DustSize()['Distribution']} \n",
        'Dust size': f"        q = {dust.get_DustSize()['q']}, a(min) = {dust.get_DustSize()['amin']} micron, a(max) = {dust.get_DustSize()['amax']} micron \n",
        'Sublimation temperature': f"        Tsub = {dust.get_Sublimation()} K \n",
        'Density Distribution': density,
        'Opacity': f"        - tau(min) = {dust.get_tau()} ; tau(max) = {dust.get_tau()} % for visual wavelength \n",
    }

def change_parameter(Path: str, change: Dict[str, str], car: Dict[str, str], nstar: int) -> None:
    """Modifies the Dusty input file based on a change dictionary."""
    file = load_file(Path)
    file = [line for line in file if ('.nk' not in line) and ('Td' not in line)]
    file = [line for line in file if all(p not in line for p in ['number of powers', 'power =', 'shell', 'Y ='])]

    if 'engelke_marengo' in change.get('Spectral', ''):
        change.pop('BB', None); change.pop('Luminosities', None)
        file = [line for line in file if 'Number of BB' not in line]
    elif 'black_body' in change.get('Spectral', ''):
        change.pop('Absorption', None)
        file = [line for line in file if 'SiO absorption depth' not in line]

    if nstar == 1 and 'Luminosities' in change:
        file = [line for line in file if 'Luminosities' not in line]
        change.pop('Luminosities')

    for param, new_line in change.items():
        idx = search_line(file, car[param])
        file[idx] = new_line
    save_file(Path, file)

def set_change(dusty, change: Dict[str, Any]) -> None:
    """Applies high-level parameter changes to the Dusty model object."""
    for key, value in change.items():
        if 'Temp' in key and key != 'Temperature':
            dusty.get_Model().get_Stars()[int(key.split('Temp')[-1])-1].set_Temperature(value)
        elif 'Lum' in key:
            dusty.get_Model().get_Stars()[int(key.split('Lum')[-1])-1].set_Luminosity(value)
        elif key == 'Opacity': dusty.get_Model().get_Dust().set_tau(value)
        elif key == 'DustSize': dusty.get_Model().get_Dust().set_DustSize(value)
        elif key == 'Sublimation': dusty.get_Model().get_Dust().set_Sublimation(value)
        elif key == 'Absorption': dusty.get_Model().set_SiOAbsorption(value)
        elif key == 'Temperature': dusty.get_Model().get_Dust().set_Temperature(value)
        elif key == 'Density': dusty.get_Model().get_Dust().set_Density(value)

def create_spectral_file(dusty, p) -> None:
    """
    Create a spectral file for the Dusty model.
    Parameters:
    dusty (Dusty): The Dusty model object.
    """

    file_name = dusty.get_Model().get_GridFile()
    bool_file,Wavelength,Flux = [], None, None
    if os.path.exists(file_name):
        for i in range(dusty.get_Model().get_NbStar()):
            if Wavelength is None:
                bool, Wavelength, Flux = get_table_interpolated(teff=p[f'Temp{i+1}'], logg=dusty.get_Model().get_Stars()[i].get_Logg(), ebv=0, gridfilename=file_name)
                bool_file.append(bool)
            else:
                bool, Wavelength, Flux_add = get_table_interpolated(teff=p[f'Temp{i+1}'], logg=dusty.get_Model().get_Stars()[i].get_Logg(), ebv=0, gridfilename=file_name)
                bool_file.append(bool)
                Flux += Flux_add
    else:
        raise FileNotFoundError(f'This file does not exist: {file_name}')
    
    if all(bool_file):
        spectral_file = dusty.get_Model().get_SpectralFile()
        with open(spectral_file, 'w') as f:
            f.write('#>')
            for i in range(len(Wavelength)):
                f.write(f'{Wavelength[i]} {Flux[i]}\n')

    else:
        spectral_file = dusty.get_Model().get_SpectralFile()
        with open(spectral_file, 'w') as f:
            f.write('#>')
            for i in range(len(Wavelength)):
                f.write(f'{Wavelength[i]} {np.nan}\n')

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



# =============================================================================
# UNIT CONVERSIONS & MATH
# =============================================================================

def list_to_dict(keys: List[str], values: List[Any]) -> Dict[str, Any]:
    """Combines two lists into a dictionary."""
    return dict(zip(keys, values))

def str_to_data(data_lines: List[str]) -> np.ndarray:
    """Converts a list of strings (space-separated) into a NumPy array."""
    return np.array([list(map(float, line.split())) for line in data_lines])

def get_column_spectrum(file_lines: List[str], index: int, index_header: int = 0) -> np.ndarray:
    """Extracts a specific column from a structured spectral file."""
    array = np.asarray(file_lines[index_header:])
    return np.array([[x for x in el.split('  ')[1:-1] if x != ''] for el in array], dtype=float).T[index]

def watt_to_jansky(Flux: np.ndarray, Wavelength: np.ndarray) -> np.ndarray:
    """Converts flux from W/m^2 to Jansky."""
    return (Flux * u.W/u.m**2).to(u.Jy, equivalencies=u.spectral_density(Wavelength * u.um)).value

def jansky_to_watt(Flux: np.ndarray, Wavelength: np.ndarray) -> np.ndarray:
    """Converts flux from Jansky to W/m^2."""
    return (Flux * u.Jy).to(u.W/u.m**2, equivalencies=u.spectral_density(Wavelength * u.um)).value

def calculate_true_radius(results: Dict[str, Any], L: float) -> float:
    """Calculates the physical radius in meters based on Dusty output and luminosity."""
    return results['r1(cm)'] * np.sqrt(L / 1e4) * 1e-2

def calculate_total_flux(F: float, r_true: float, distance: float) -> float:
    """Calculates total flux at a given distance."""
    return F * (r_true**2) / (distance**2)

def interpolate_spline(Wavelength: np.ndarray, Flux: np.ndarray, order: int = 3):
    """Creates a spline interpolation for flux data."""
    return make_interp_spline(Wavelength, Flux, k=order)

def interpol(x_eval, x_model, y_model):
    """Interpolates model values onto an evaluation grid."""
    try:
        return interpolate_spline(np.asarray(x_model).flatten(), np.asarray(y_model).flatten())(np.asarray(x_eval).flatten())
    except:
        return interpolate_spline(x_model, y_model)(x_eval)
    
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



# =============================================================================
# MCMC MODEL WRAPPERS (PyMCMC Compatible)
# =============================================================================

def model(theta: np.ndarray, data: Any) -> np.ndarray:
    """
    Main model wrapper for PyMCMC. Updates Dusty, runs it, and returns interpolated flux.
    """
    dusty_obj, data_mod, fit_obj, logfile, jansky, lock = data.user_defined_object[0]
    p = dict(zip(fit_obj.get_Param().keys(), theta))
    
    # Update Dust Properties
    dust_size = dusty_obj.get_Model().get_Dust().get_DustSize()
    density = dusty_obj.get_Model().get_Dust().get_Density()
    if 'amin' in p: dust_size['amin'] = np.round(10**p['amin'], 3)
    if 'amax' in p: dust_size['amax'] = np.round(10**p['amax'], 3)
    if 'q' in p: dust_size['q'] = p['q']
    if 'shell' in p: density['shell'] = p['shell']

    # Validity checks
    if dust_size['amin'] > dust_size['amax'] or density['shell'] < 0:
        return np.full_like(data.xdata, np.nan)

    # Apply changes and run Dusty
    change = {k: v for k, v in p.items() if k not in ['amin', 'amax', 'q', 'shell']}
    change.update({'DustSize': dust_size, 'Density': density})
    set_change(dusty_obj, change)
    dusty_obj.change_parameter()

    if dusty_obj.get_Model().get_Spectral() in ['FILE_LAMBDA_F_LAMBDA', 'FILE_F_LAMBDA']:
        create_spectral_file(dusty_obj, p)

    # Use lock if running in parallel mode
    if lock is not None:
        with lock:
            dusty_obj.lunch_dusty(verbose=0, logfile=logfile)
            dusty_obj.make_SED(luminosity=p['Lest'], Jansky=jansky)
    else:
        dusty_obj.lunch_dusty(verbose=0, logfile=logfile)
        dusty_obj.make_SED(luminosity=p['Lest'], Jansky=jansky)

    # Bandpass integration or interpolation
    if data_mod.get_table() is not None:
        bandpass_dict = data_mod.get_common_filters()
        y_model, _ = dusty_obj.get_SED().integrate_bandpass(bandpass=bandpass_dict)
        # Ensure correct order
        pivot_waves = [get_central_wavelegnth(get_bandpass(f)) for f in bandpass_dict.values()]
        y_model = y_model[np.argsort(pivot_waves)]
    else:
        y_model = interpol(data.xdata, dusty_obj.get_SED().get_Wavelength(), dusty_obj.get_SED().get_Flux())
        
    return y_model

def prediction_model(theta, context, x_eval):
    """Model wrapper specifically for generating prediction intervals across a grid."""
    dusty_obj, data_mod, fit_obj, logfile, jansky = context
    p = dict(zip(fit_obj.get_Param().keys(), theta))
    # ... (similar logic to model function, updates dusty_obj and returns interpolated flux)
    return interpol(x_eval, dusty_obj.get_SED().get_Wavelength(), dusty_obj.get_SED().get_Flux())

def chi2(theta: np.ndarray, data: Any) -> float:
    """Calculates Chi-squared between observations and the current model iteration."""
    try:
        x, y = data.xdata[0], data.ydata[0]
    except AttributeError:
        x, y = data
    y_model = interpol(theta, x, y).reshape(y.shape)
    return np.nansum((y_model - y)**2)

# =============================================================================
# ASTRONOMY & FILTERS
# =============================================================================

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

def get_bandpass(bandpass_name: str) -> SpectralElement:
    """Loads a bandpass filter from a FITS file."""
    path = glob.glob(os.path.join(os.path.dirname(__file__), 'filter', 'comp', 'nonhst', f"{bandpass_name}.fits"))
    if not path: raise FileNotFoundError(f"Filter {bandpass_name} not found.")
    return SpectralElement.from_file(path[0])

def intergrate_bandpass(wavelength: np.ndarray, flux: np.ndarray, bandpass: SpectralElement, n_int: int = 10000):
    """Integrates a spectrum through a bandpass filter using custom Simpson's rule."""
    w_int = np.linspace(np.min(bandpass.waveset.value), np.max(bandpass.waveset.value), n_int) * u.AA
    f_interp = np.interp(w_int.value / 10000, wavelength, flux)
    thru = bandpass(w_int) * f_interp
    num = simpson(w_int.value / 10000, thru)
    den = simpson(w_int.value / 10000, bandpass(w_int))
    return num / den if den != 0 else 0, 0 # Error logic omitted for brevity

def get_central_wavelegnth(bandpass: SpectralElement) -> float:
    """Returns the pivot wavelength of a bandpass."""
    return bandpass.pivot().value

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

def aggregate_table(table, column: str = 'sed_filter', fct: Callable = np.mean):
    """
    Agrège les données d'une table en fonction d'une fonction donnée.

    Paramètres:
    table (Table): La table à agréger.
    function (function): La fonction d'agrégation à appliquer.

    Retourne:
    Table: La table agrégée.
    """

    units = [table[col].unit for col in ['sed_filter', 'sed_freq', 'sed_flux', 'sed_eflux', '_tabname'] if col in table.colnames and hasattr(table[col], 'unit')]

    grouped = table.group_by('sed_filter')
    keys = grouped.groups.keys['sed_filter']
    new_rows = []
    for i, group in enumerate(grouped.groups):
        sed_freq_mean = fct(group['sed_freq'])
        sed_flux_mean = fct(group['sed_flux'])
        sed_eflux_mean = fct(group['sed_eflux'])
        tabname_list = list(group['_tabname'])
        new_rows.append((str(keys[i]), sed_freq_mean, sed_flux_mean, sed_eflux_mean, tabname_list))

    # Create the new table
    return Table(rows=new_rows, names=['sed_filter','sed_freq', 'sed_flux', 'sed_eflux', 'tabname'], units=units)


#==============================================================================
# Bandpass creation for custom spectra
#==============================================================================

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

if __name__ == "__main__":
    pass