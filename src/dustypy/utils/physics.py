"""
Physics and Mathematical Utilities module for DustyPY.
Includes integration routines, extinction laws, grid interpolation, 
and wavelength grid management.
"""

import os
import ctypes
from typing import Dict, List, Optional, Tuple, Any

import urllib.request
from pathlib import Path

import numpy as np
from scipy import interpolate
from scipy.interpolate import LinearNDInterpolator
from astropy.io import fits
import astropy.units as u


# =============================================================================
# C LIBRARY LOADING (Simpson's Rule)
# =============================================================================

class SimpsonResult(ctypes.Structure):
    """
    C-compatible structure to receive results from the Simpson shared library.
    
    Attributes:
        value (float): The calculated integral value.
        error (float): Estimated numerical error.
    """
    _fields_ = [("value", ctypes.c_double),
                ("error", ctypes.c_double)]


def _load_simpson_library() -> Optional[ctypes.CDLL]:
    """
    Dynamically loads the compiled Simpson shared library from the package.

    Searches the 'libs' directory for a shared object file starting with 
    'simpson'.

    Returns:
        Optional[ctypes.CDLL]: The loaded C library or None if not found.
    """
    lib_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "libs")
    if not os.path.exists(lib_dir):
        return None

    for file in os.listdir(lib_dir):
        if file.startswith("simpson") and file.endswith(".so"):
            try:
                lib = ctypes.CDLL(os.path.join(lib_dir, file))
                lib.simpson.argtypes = [
                    ctypes.POINTER(ctypes.c_double), 
                    ctypes.POINTER(ctypes.c_double), 
                    ctypes.c_int
                ]
                lib.simpson.restype = SimpsonResult
                return lib
            except (OSError, AttributeError):
                continue
    return None


_SIMPSON_LIB = _load_simpson_library()


def simpson_integrate(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """
    Integrates a function using Simpson's rule via an optimized C library.

    Falls back to scipy.integrate.simpson if the C library is unavailable.

    Args:
        x (np.ndarray): Wavelength/independent variable grid.
        y (np.ndarray): Flux/dependent variable values.

    Returns:
        Tuple[float, float]: A tuple containing (Integral Value, Estimated Error).
            Error is 0.0 if the Scipy fallback is used.
    """
    if _SIMPSON_LIB is None:
        from scipy.integrate import simpson
        return float(simpson(y=y, x=x)), 0.0

    # Simpson's rule requires an odd number of points
    n = len(x) if len(x) % 2 == 1 else len(x) - 1
    
    x_c = np.ascontiguousarray(x[:n], dtype=np.double).ctypes.data_as(
        ctypes.POINTER(ctypes.c_double)
    )
    y_c = np.ascontiguousarray(y[:n], dtype=np.double).ctypes.data_as(
        ctypes.POINTER(ctypes.c_double)
    )

    res = _SIMPSON_LIB.simpson(x_c, y_c, n)
    return res.value, res.error


# =============================================================================
# PHYSICS UTILS
# =============================================================================

def log_space(start: float, stop: float, num: int) -> np.ndarray:
    """
    Creates a logarithmically spaced array.

    Args:
        start (float): Start value.
        stop (float): Stop value.
        num (int): Number of points.

    Returns:
        np.ndarray: The generated log-grid.
    """
    return np.logspace(np.log10(start), np.log10(stop), num)


def validate_wavelength_spacing(wavelengths: np.ndarray) -> bool:
    """
    Validates that wavelength spacing does not exceed Dusty's 50% stability limit.

    Args:
        wavelengths (np.ndarray): The wavelength grid to check.

    Returns:
        bool: True if validation passes.

    Raises:
        ValueError: If spacing between two consecutive points is >= 50%.
    """
    ratios = wavelengths[1:] / wavelengths[:-1]
    if not np.all(ratios < 1.5):
        idx = np.where(ratios >= 1.5)[0][0]
        raise ValueError(
            f"Wavelength spacing too large (>50%) at index {idx}: "
            f"{wavelengths[idx]} -> {wavelengths[idx+1]}. "
            "Dusty requires smoother transitions for convergence."
        )
    return True


def apply_extinction(
    wavelength_um: np.ndarray, 
    flux: np.ndarray, 
    ebv: float, 
    rv: float = 3.1, 
    lmc2: bool = False, 
    avglmc: bool = False
) -> np.ndarray:
    """
    Applies interstellar extinction correction based on Fitzpatrick (1999).

    Args:
        wavelength_um (np.ndarray): Wavelength array in microns.
        flux (np.ndarray): Input flux array.
        ebv (float): E(B-V) color excess. Positive values redden the spectrum, 
            negative values de-redden it.
        rv (float): Total-to-selective extinction ratio. Defaults to 3.1.
        lmc2 (bool): Use Large Magellanic Cloud (LMC2) parameters.
        avglmc (bool): Use average LMC parameters.

    Returns:
        np.ndarray: The extinction-corrected flux.
    """
    # Early exit if extinction is negligible
    if np.abs(ebv) < 1e-10:
        return flux

    x = 1. / wavelength_um  # Inverse microns
    curve = np.zeros_like(x)

    # Standard Fitzpatrick (1999) parameters
    params = {
        "x0": 4.596, "gamma": 0.99, "c3": 3.23, "c4": 0.41,
        "c2": -0.824 + 4.717 / rv, 
        "c1": 2.030 - 3.007 * (-0.824 + 4.717 / rv)
    }

    if lmc2:
        params.update({
            "x0": 4.626, "gamma": 1.05, "c3": 1.92, 
            "c4": 0.42, "c2": 1.31, "c1": -2.16
        })
    elif avglmc:
        params.update({
            "x0": 4.596, "gamma": 0.91, "c3": 2.73, 
            "c4": 0.64, "c2": 1.11, "c1": -1.28
        })

    # UV portion
    x_cut_uv = 10000.0 / 2700.0
    i_uv = np.where(x >= x_cut_uv)[0]
    
    # Calculate UV curve for transition spline points
    x_uv_ref = np.array([10000.0 / 2700.0, 10000.0 / 2600.0])
    x_uv_calc = np.concatenate((x_uv_ref, x[i_uv])) if len(i_uv) > 0 else x_uv_ref
    
    # UV polynomial + Lorentzian bump
    y_uv = (params["c1"] + params["c2"] * x_uv_calc + 
            params["c3"] * x_uv_calc**2 / 
            ((x_uv_calc**2 - params["x0"]**2)**2 + (x_uv_calc * params["gamma"])**2))
    
    # UV far-UV curvature
    y_uv += params["c4"] * (
        0.5392 * (np.maximum(x_uv_calc, 5.9) - 5.9)**2 + 
        0.05644 * (np.maximum(x_uv_calc, 5.9) - 5.9)**3
    ) + rv
    
    if len(i_uv) > 0:
        curve[i_uv] = y_uv[2:]

    # Optical/IR portion via spline interpolation
    x_spline = np.concatenate(
        ([0], 10000.0 / np.array([26500.0, 12200.0, 6000.0, 5470.0, 4670.0, 4110.0]))
    )
    y_spline = np.concatenate((
        np.array([0.0, 0.26469, 0.82925]) * rv / 3.1, 
        [
            np.polyval([-0.422809, 1.0027, 0.000213572][::-1], rv),
            np.polyval([-0.051354, 1.00216, -0.000073578][::-1], rv),
            np.polyval([0.700127, 1.00184, -0.00003326][::-1], rv),
            np.polyval([1.19456, 1.01707, -0.00546959, 0.000797809, -0.0000445636][::-1], rv)
        ]
    ))
    
    i_op_ir = np.where(x < x_cut_uv)[0]
    if len(i_op_ir) > 0:
        tck = interpolate.splrep(
            np.concatenate((x_spline, x_uv_ref)), 
            np.concatenate((y_spline, y_uv[:2])), 
            s=0
        )
        curve[i_op_ir] = interpolate.splev(x[i_op_ir], tck)

    # Numerical Safety: Clip the exponent to prevent overflow in Chi2 calculations
    exponent = 0.4 * curve * ebv
    exponent = np.clip(exponent, -100, 100)
    
    return flux * 10.**exponent


# =============================================================================
# FILTER UTILS
# =============================================================================

def get_common_filters(vizier_filters: List[str], local_filenames: List[str]) -> Dict[str, str]:
    """
    Matches VizieR SED filter names with local package bandpass filenames.

    Example:
        vizier_filters: ['2MASS:J', 'WISE:W1']
        local_filenames: ['2mass_j.fits', 'wise_w1.fits']
        Returns: {'2MASS:J': '2mass_j', 'WISE:W1': 'wise_w1'}

    Args:
        vizier_filters (List[str]): List of filter strings from VizieR tables.
        local_filenames (List[str]): List of available .fits files in the package.

    Returns:
        Dict[str, str]: Mapping from VizieR names to local filenames.
    """
    available_bp = {
        os.path.splitext(f)[0].lower(): os.path.splitext(f)[0] 
        for f in local_filenames
    }
    
    common_filters = {}
    for viz_filter in vizier_filters:
        normalized_f = viz_filter.lower().replace(':', '_')
        
        for bp_key, original_name in available_bp.items():
            if bp_key in normalized_f:
                common_filters[viz_filter] = original_name
                break
                
    return common_filters


# =============================================================================
# ATMOSPHERE UTILS
# =============================================================================

_ATMOSPHERE_GRID_CACHE: Dict[str, Dict[str, Any]] = {}


def get_interpolated_atmosphere(
    grid_path: str, 
    teff: float, 
    logg: float
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """
    Loads a FITS atmosphere grid and interpolates flux for specific parameters.

    Results are cached in memory to speed up repeated calls during fitting.

    Expected FITS structure:
    - HDU 0: Header 'WAVUNIT' (unit name).
    - HDU 1-N: Headers 'TEFF', 'LOGG' (parameters).
    - HDU 1-N: Table with 'flux' column.
    - HDU 1: Table with 'wavelength' column (used for all HDUs).

    Args:
        grid_path (str): Path to the FITS atmosphere grid.
        teff (float): Effective temperature in Kelvin.
        logg (float): Surface gravity (log10).

    Returns:
        Tuple[np.ndarray, Optional[np.ndarray]]: A tuple containing 
            (Wavelength array, Flux array). Flux is None if parameters are 
            out of grid bounds.
    """
    global _ATMOSPHERE_GRID_CACHE
    
    if grid_path not in _ATMOSPHERE_GRID_CACHE:
        if not os.path.exists(grid_path):
            try:
                download_atmosphere_grid(
                    destination_path=grid_path,
                    zenodo_id="",  # Example Zenodo ID for MARCS grid
                    filename=os.path.basename(grid_path)
                )
            except Exception:
                raise FileNotFoundError(
                    f"Atmosphere grid file not found at {grid_path} "
                    "and automatic download failed."
                )
            
        with fits.open(grid_path) as hdu:
            # 1. Validate basic structure
            if len(hdu) < 2:
                raise ValueError(f"FITS file {grid_path} must have at least one data extension.")
            
            wave_unit_str = hdu[0].header.get('WAVUNIT', 'micron')
            wave_unit = u.Unit(wave_unit_str)
            
            teffs, loggs, fluxes = [], [], []
            
            # 2. Extract Wavelength from the first extension
            if 'Wavelength' not in hdu[1].columns.names:
                raise KeyError(f"Missing 'Wavelength' column in HDU 1 of {grid_path}")
            wavelength = (hdu[1].data['Wavelength'] * wave_unit).to(u.um).value
            
            # 3. Extract Grid Points
            for i in range(1, len(hdu)):
                header = hdu[i].header
                # Check for mandatory headers
                if 'TEFF' not in header or 'LOGG' not in header:
                    print(f"Warning: Skipping HDU {i}, missing TEFF or LOGG.")
                    continue
                
                teffs.append(header['TEFF'])
                loggs.append(header['LOGG'])
                fluxes.append(hdu[i].data['Flux'])
            
            if not teffs:
                raise ValueError("No valid TEFF/LOGG data found in the grid.")
                
            # Setup 2D Linear Interpolator
            points = np.array(list(zip(teffs, loggs)))
            flux_grid = np.array(fluxes)
            interpolator = LinearNDInterpolator(points, flux_grid)
            
            _ATMOSPHERE_GRID_CACHE[grid_path] = {
                'wavelength': wavelength,
                'interpolator': interpolator,
                'teff_range': (float(min(teffs)), float(max(teffs))),
                'logg_range': (float(min(loggs)), float(max(loggs)))
            }

    grid_data = _ATMOSPHERE_GRID_CACHE[grid_path]
    
    # 4. Bound checking
    t_min, t_max = grid_data['teff_range']
    g_min, g_max = grid_data['logg_range']
    
    if not (t_min <= teff <= t_max) or not (g_min <= logg <= g_max):
        # We return None for flux if outside bounds to allow the Runner 
        # to handle the fallback (e.g., using a Blackbody instead).
        return grid_data['wavelength'], None
        
    flux = grid_data['interpolator'](teff, logg)
    
    if np.isnan(flux).any():
        return grid_data['wavelength'], None
        
    return grid_data['wavelength'], flux

def download_atmosphere_grid(destination_path: str, zenodo_id: str, filename: str):
    """
    Télécharge une grille d'atmosphère depuis Zenodo si elle est absente.
    
    Args:
        destination_path (str): Chemin local où enregistrer le fichier.
        zenodo_id (str): L'identifiant de votre dépôt Zenodo (ex: '1234567').
        filename (str): Le nom du fichier sur Zenodo (ex: 'marcs_z0.00.fits').
    """
    path = Path(destination_path)
    if path.exists():
        return # Le fichier est déjà là

    # URL de téléchargement direct Zenodo
    url = f"https://zenodo.org/record/{zenodo_id}/files/{filename}?download=1"
    
    print(f"Downloading {filename} from Zenodo (~800MB)...")
    
    # Création du dossier parent si nécessaire
    path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        urllib.request.urlretrieve(url, destination_path)
        print(f"Download complete: {destination_path}")
    except Exception as e:
        print(f"Error downloading from Zenodo: {e}")
        raise


def get_grid_bounds(grid_path: str) -> Tuple[Optional[Tuple[float, float]], Optional[Tuple[float, float]]]:
    """
    Retrieves the Teff and Logg boundaries of a cached atmosphere grid.

    Args:
        grid_path (str): Path to the grid file.

    Returns:
        Tuple[Optional[Tuple[float, float]], Optional[Tuple[float, float]]]: 
            (Teff min/max, Logg min/max). Returns (None, None) if grid not cached.
    """
    if grid_path in _ATMOSPHERE_GRID_CACHE:
        data = _ATMOSPHERE_GRID_CACHE[grid_path]
        return data['teff_range'], data['logg_range']
    return None, None