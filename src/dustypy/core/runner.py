"""
Dusty execution runner module.
Handles the orchestration of the Dusty binary execution, environment isolation, 
and output parsing.
"""

import os
import subprocess
from typing import Dict, Optional, Tuple

import numpy as np

from ..core.model import Model
from ..utils.physics import get_interpolated_atmosphere
from ..io.writer import Writer
from .result import Result


class Runner:
    """
    Manages the execution environment and running of the Dusty binary.

    This class handles the creation of isolated execution directories to 
    prevent file collisions during parallel runs, manages symlinks to 
    Dusty's required data libraries, and parses the resulting output files.

    Attributes:
        binary_path (str): Absolute path to the Dusty executable.
        bin_dir (str): Directory containing the Dusty binary and data folders.
    """

    def __init__(self, binary_path: str):
        """
        Initializes the Runner.

        Args:
            binary_path (str): System path to the compiled Dusty binary.
        """
        self.binary_path = os.path.abspath(binary_path)
        self.bin_dir = os.path.dirname(self.binary_path)

    def run(self, model: Model, run_dir: str = "tmp_run") -> Result:
        """
        Executes Dusty in an isolated environment and manages model state.

        This method prepares the execution folder, generates the necessary 
        input files, runs the binary via subprocess, and parses the 
        results into a Result object. It ensures that the model state is 
        restored after execution, which is critical for MCMC iterations.

        Args:
            model (Model): The physical model configuration to run.
            run_dir (str): Relative or absolute path for the temporary 
                execution directory. Defaults to "tmp_run".

        Returns:
            Result: An object containing the parsed SED and scalar outputs.

        Raises:
            subprocess.CalledProcessError: If the Dusty binary fails during execution.
            FileNotFoundError: If required Dusty data directories are missing.
        """
        run_path = os.path.abspath(run_dir)
        os.makedirs(run_path, exist_ok=True)

        # 1. Environment Isolation (Symlinking required folders)
        # Link the optical properties library (Lib_nk)
        lib_src = os.path.join(self.bin_dir, 'Lib_nk')
        lib_dst = os.path.join(run_path, 'Lib_nk')
        if os.path.exists(lib_src) and not os.path.exists(lib_dst):
            os.symlink(lib_src, lib_dst)

        # Create a local data directory and link standard data files
        data_dst = os.path.join(run_path, 'data')
        os.makedirs(data_dst, exist_ok=True)
        
        data_src = os.path.join(self.bin_dir, 'data')
        if os.path.exists(data_src):
            for item in os.listdir(data_src):
                # Skip the global wavelength grid if a custom one is provided
                if item == "lambda_grid.dat" and model.wavelengths is not None:
                    continue
                s = os.path.join(data_src, item)
                d = os.path.join(data_dst, item)
                if not os.path.exists(d):
                    os.symlink(s, d)

        # 2. Spectral State Management
        # Save original state to ensure restoration after the run
        original_shape = model.spectral_shape
        original_file = model.spectral_file

        try:
            # Handle Atmosphere Grid mode
            if original_shape == "atmosphere":
                # Generate the composite spectral file from the FITS grid
                generated_file = self._prepare_spectral_file(model, run_path)
                # Temporarily switch Dusty to file input mode
                model.spectral_shape = "file_lambda_f_lambda"
                model.spectral_file = generated_file

            # 3. Input Generation
            writer = Writer(model)
            writer.write_input_file(run_path)
            
            # Write custom wavelength grid if defined
            if model.wavelengths is not None:
                original_grid_template = os.path.join(
                    self.bin_dir, "data", "lambda_grid.dat"
                )
                writer.write_wavelength_file(run_path, original_grid_template)

            # 4. Create the Dusty pointer file (.mas)
            with open(os.path.join(run_path, "model.mas"), 'w') as f:
                f.write(f"{model.name}.inp\n")

            # 5. Execute Dusty Binary
            subprocess.run(
                [self.binary_path, "model.mas"], 
                cwd=run_path, 
                check=True, 
                capture_output=True
            )

            # 6. Result Parsing
            stb_path = os.path.join(run_path, f"{model.name}.stb")
            out_path = os.path.join(run_path, f"{model.name}.out")
            
            # Load SED (Wavelength vs lambda*F_lambda)
            spectrum_data = np.loadtxt(stb_path, skiprows=63)
            # Parse scalar parameters from the .out file
            scalars = self._parse_out_file(out_path)

            # 7. Create Result Object and compute synthetic photometry
            result = Result(
                wavelength=spectrum_data[:, 0],
                flux_norm=spectrum_data[:, 1],
                scalars=scalars,
                model=model
            )
            result.compute_all_photometry(model.dataset)
            
            return result

        except subprocess.CalledProcessError as e:
            # Log binary error output for debugging
            print(f"Dusty Binary Execution Error:\n{e.stderr.decode()}")
            raise
        
        finally:
            # RESTORE MODEL STATE
            # Critical for MCMC to allow the next iteration to enter 
            # the 'atmosphere' logic block again.
            model.spectral_shape = original_shape
            model.spectral_file = original_file

    def _parse_out_file(self, file_path: str) -> Dict[str, float]:
        """
        Parses the scalar results table from the Dusty .out file.

        Args:
            file_path (str): Path to the .out file.

        Returns:
            Dict[str, float]: A dictionary mapping parameter names to 
                their calculated values.
        """
        with open(file_path, 'r') as f:
            lines = f.readlines()

        header_idx = -1
        for i, line in enumerate(lines):
            if 'tau0' in line:
                header_idx = i
                break
        
        if header_idx == -1:
            return {}

        # Extract keys from the header line (stripping comment characters)
        keys = lines[header_idx].replace('###', '').split()
        # The values are typically 3 lines below the header
        values_line = lines[header_idx + 3].split()
        # Skip the first element which is usually the model index
        values = [float(v) for v in values_line[1:]]
        
        return dict(zip(keys, values))
    
    def _prepare_spectral_file(self, model: Model, run_path: str) -> str:
        """
        Generates a composite spectrum file from an atmosphere grid.

        Loops through all stars in the model, interpolates their spectra 
        from the FITS grid, and sums them into a single spectral file 
        readable by Dusty.

        Args:
            model (Model): The model containing star parameters and grid path.
            run_path (str): The current execution directory.

        Returns:
            str: The filename of the generated spectral data.

        Raises:
            ValueError: If the grid path is missing or parameters are 
                outside the grid range.
        """
        if model.grid_path is None:
            raise ValueError(
                "model.grid_path must be provided to use atmosphere models."
            )

        composite_flux = None
        common_wavelength = None

        for star in model.stars:
            # Perform 2D interpolation (Teff, logg) via physics utility
            wl, flux = get_interpolated_atmosphere(
                model.grid_path, 
                star['temperature'], 
                star['logg']
            )
            
            if flux is None:
                raise ValueError(
                    f"Star (T={star['temperature']}, logg={star['logg']}) "
                    f"is outside the atmosphere grid range."
                )

            # Sum the fluxes weighted by their absolute luminosities
            # Dusty will handle final total normalization.
            star_flux = flux * star['luminosity']
            
            if composite_flux is None:
                common_wavelength = wl
                composite_flux = star_flux
            else:
                composite_flux += star_flux

        spec_filename = "atmosphere_composite.dat"
        spec_full_path = os.path.join(run_path, spec_filename)
        
        # Write file in Dusty's '#> lambda flux' format
        with open(spec_full_path, 'w') as f:
            f.write("#> Composite atmosphere spectrum generated by DustyPY\n")
            for w, f_val in zip(common_wavelength, composite_flux):
                f.write(f"{w:.6e} {f_val:.6e}\n")
        
        return spec_filename