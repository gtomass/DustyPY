# DustyPY

**DustyPY** is a modern, high-performance Python wrapper for the Dusty radiative transfer code. It is designed to bridge the gap between legacy Fortran power and modern data science workflows, enabling automated parameter estimation of stellar envelopes through Markov Chain Monte Carlo (MCMC) sampling.



## Main Features

* **Pythonic Interface**: Full object-oriented control over Dusty models (Stars, Dust, Geometry).

* **MCMC Integration**: Seamless coupling with `PyMCMC` for automated fitting of temperatures, opacities, and distances.

* **Atmosphere Grids**: Support for MARCS and other FITS atmosphere grids with high-speed 2D interpolation.

* **Observational Data**: Automatic SED retrieval from **VizieR** and intelligent photometric aggregation.

* **Synthetic Photometry**: Fast integration of model spectra over instrument bandpasses using a cached filter system.

* **Performance**: C-optimized integration routines (Simpson's rule) with **OpenMP** multi-threading support.


## Installation

**Prerequisites**

1. **Dusty Binary**: You must have the [`dusty`](https://github.com/ivezic/dusty) executable compiled on your system.

2. **C Compiler**: A compiler supporting OpenMP (gcc, clang, or msvc) is required to build the performance extensions.

**Install from Source**

Clone the repository and install it in editable mode:

```bash
git clone https://github.com/gtomassini/DustyPY.git
cd DustyPY
pip install -e .
```

This will automatically fetch dependencies like `numpy`, `astropy`, and the latest version of `PyMCMC` from GitHub.

## Quick Start

**1. Define a Model**

```python
from dustypy import Model, Runner

# Setup a Betelgeuse-like model
model = Model(name="Betelgeuse", distance=222.0)
model.add_star(temperature=3600, luminosity=1.2e5, logg=-0.5)
model.set_dust(opacity=0.1, temperature=400)
```

**2. Run dusty**
```python
runner = Runner(binary_path="/path/to/dusty")
result = runner.run(model)

# Access physical results
print(result) # Beautiful summary of radii, mass-loss, and photometry
result.plot() # Instant SED visualization
```

**3. MCMC Fitting**
```python
from dustypy import Fitter, Dataset

# Get real data from VizieR
ds = Dataset.from_vizier("Betelgeuse")

fitter = Fitter(model, runner, ds, n_workers=4)
fitter.add_star_parameter(0, 'temp', initial=3600, min_val=3000, max_val=4000)
fitter.add_parameter('opacity', initial=0.1, min_val=0.0, max_val=1.0)

analyzer = fitter.run(n_iterations=1000)
analyzer.plot_traces()
analyzer.plot_corner()
fitter.best_model.plot()
```

## Project Structure

* `src/dustypy/core/`: Core logic (Model, Runner, Result, Fitter, Likelihood).

* `src/dustypy/utils/`: Physics utilities and C-wrappers.

* `src/dustypy/io/`: Input/Output management and templates.

* `src/dustypy/data/`: Instrumental filters and extinction laws.

## Note

As of today, there is a bug in dusty V4 that has not yet been corrected. This implies that the is not taken into account. Refer to the [raised issue](https://github.com/ivezic/dusty/issues/11) for solution.

## Data Specifications

### Atmosphere Grids (FITS)
To ensure compatibility with the interpolation engine, FITS grids must follow this structure:
* **HDU 0**: Primary header. Optional `WAVUNIT` keyword (e.g., 'Angstrom', 'micron'). Defaults to 'micron'.
* **HDU 1**: Contains the reference `wavelength` column and the first `flux` column.
* **HDU 1 to N**: Each extension must contain:
    * **Header**: `TEFF` (float) and `LOGG` (float) keywords.
    * **Data**: A table with at least a `flux` column.
* **Note**: The code assumes all extensions share the same wavelength grid as HDU 1.

You can find Marcs grid in the good format on [Zenodo](https://zenodo.org/uploads/18412489) or download them directly in the code using:

```python
from dustypy import utils

PATH = 'PATH/TO/LOCATION/YOU/WANT'
file = 'marcs_z0.00.fits'

utils.download_atmosphere_grid(destination_path=PATH, filename=file)

```

## License

DustyPY is licensed under the MIT License. See the [LICENSE](https://github.com/Radron74/DustyPY/blob/main/LICENSE) file for more details.