"""
Dusty Input Writer module.
Handles the generation of .inp files and custom wavelength grids by mapping 
Python Model objects to the specific text format required by Dusty.
"""

import os
from typing import Dict, Optional, Any
from importlib import resources

from ..core.model import Model


class Writer:
    """
    Handles the generation of Dusty input files via template replacement.

    This class reads a base template file containing placeholders and populates 
    them with physical parameters from a Model instance. It also manages the 
    creation of custom wavelength grid files.

    Attributes:
        model (Model): The model instance providing physical parameters.
    """

    def __init__(self, model: Model):
        """
        Initializes the Writer with a specific model.

        Args:
            model (Model): The model configuration to write.
        """
        self.model = model
        # Use importlib.resources to find the template within the package
        self._template_ref = resources.files('dustypy').joinpath('io/templates/Mod.inp')

    def write_input_file(self, destination_path: str, filename: Optional[str] = None) -> str:
        """
        Generates the .inp file in the specified directory.

        Args:
            destination_path (str): Directory where the file will be saved.
            filename (Optional[str]): Name of the output file. 
                Defaults to {model.name}.inp.

        Returns:
            str: The absolute path to the generated .inp file.

        Raises:
            FileNotFoundError: If the Mod.inp template is missing from the package.
        """
        if filename is None:
            filename = f"{self.model.name}.inp"
            
        if not self._template_ref.is_file():
            raise FileNotFoundError(f"Dusty template not found at {self._template_ref}")

        # Read the template content
        with open(self._template_ref, 'r', encoding='utf-8') as f:
            content = f.read()

        # Generate the replacement dictionary based on model state
        replacements = self._generate_replacements()
        
        # Replace placeholders in the template
        for key, value in replacements.items():
            content = content.replace(key, str(value))

        # Ensure the output directory exists
        os.makedirs(destination_path, exist_ok=True)
        final_path = os.path.join(destination_path, filename)
        
        with open(final_path, 'w', encoding='utf-8') as f:
            f.write(content)
        
        return final_path
    
    def write_wavelength_file(self, run_dir: str, template_grid_path: str):
        """
        Generates a custom lambda_grid.dat file for the simulation.

        If the model has a custom wavelength grid defined, this method 
        overwrites the default grid link in the run directory with a 
        formatted text file.

        Args:
            run_dir (str): The temporary directory for the Dusty run.
            template_grid_path (str): Path to the original lambda_grid.dat 
                to preserve header comments.
        """
        if self.model.wavelengths is None:
            return  # Use default grid if none is specified

        target_path = os.path.join(run_dir, "data", "lambda_grid.dat")
        
        # Remove symlink if it exists to replace it with a real file
        if os.path.islink(target_path):
            os.unlink(target_path)

        # 1. Read original header to preserve comments and format
        header = []
        if os.path.exists(template_grid_path):
            with open(template_grid_path, 'r', encoding='utf-8') as f:
                for line in f:
                    if '# nL =' in line:
                        header.append(f'# nL = {len(self.model.wavelengths)}\n')
                    elif line.startswith('#'):
                        header.append(line)
                    else:
                        break

        # 2. Write the new custom grid
        os.makedirs(os.path.dirname(target_path), exist_ok=True)
        with open(target_path, 'w', encoding='utf-8') as f:
            f.writelines(header)
            for w in self.model.wavelengths:
                f.write(f"{w:.6e}\n")

    def _generate_replacements(self) -> Dict[str, str]:
        """
        Constructs the logic strings for Dusty's input parameters.

        This method handles the complex formatting requirements of Dusty, 
        specifically managing spectral modes (BlackBody vs File) and 
        dust/geometry parameters.

        Returns:
            Dict[str, str]: A dictionary mapping template placeholders 
                to formatted strings.
        """
        model = self.model
        dust = model.dust
        geometry = model.geometry

        # Initialize placeholders as empty strings
        spectral_logic = ""
        bb_logic = ""
        temp_logic = ""
        lum_logic = ""
        abs_logic = ""

        spectral_shape = model.spectral_shape
        
        # --- 1. Spectral Mode Handling ---
        if spectral_shape == 'black_body':
            spectral_logic = "      Spectral shape = black_body"
            bb_logic = f"        Number of BB = {len(model.stars)}"
            temp_logic = f"        Temperature = {', '.join(str(s['temperature']) for s in model.stars)} K"
            
            # Multi-BB models require explicit relative luminosities
            if len(model.stars) > 1:
                total_l = model.total_luminosity
                lums_rel = [f"{s['luminosity']/total_l:.4f}" for s in model.stars]
                lum_logic = f"        Luminosities = {', '.join(lums_rel)}"

        elif spectral_shape == 'engelke_marengo':
            spectral_logic = "      Spectral shape = engelke_marengo"
            # Engelke mode assumes a single temperature source
            temp_logic = f"        Temperature = {model.stars[0]['temperature']} K"
            abs_logic = f"        SiO absorption depth = {dust.get('sio_absorption', 10.0)} percents"

        else:
            # File-based spectral modes (e.g., file_lambda_f_lambda)
            # In these modes, Dusty ignores BB-specific parameter lines
            spectral_file = model.spectral_file if model.spectral_file else "source.nk"
            # Double indentation (\t\t) is required by Dusty to parse file paths correctly
            spectral_logic = f"      Spectral shape = {spectral_shape} \n \t\t{spectral_file}"

        # --- 2. Dust Composition and Abundances ---
        comp = dust.get('composition', {})
        comp_files = "\n        ".join(f"data/Lib_nk/{name}.nk" for name in comp.keys())
        comp_logic = (
            f"        Number of additional components = {len(comp)}, "
            f"properties listed in files \n         {comp_files}"
        )
        
        abund_str = ", ".join(str(v) for v in comp.values())
        abund_logic = f"        Abundances for these components = {abund_str}"

        # --- 3. Grain Size Distribution ---
        sd = dust.get('size_distribution', {
            'type': 'MODIFIED_MRN', 'q': 3.5, 'amin': 0.005, 'amax': 0.25
        })
        sd_logic = f"        SIZE DISTRIBUTION = {sd['type']}"
        sd_params = (f"        q = {sd['q']}, a(min) = {sd['amin']} micron, "
                     f"a(max) = {sd['amax']} micron")

        # --- 4. Density Profile and Envelope Geometry ---
        if geometry['type'].upper() == 'POWD':
            density_logic = (
                f"         density type = POWD\n"
                f"        \t number of powers = {geometry.get('n_powers', 1)}\n"
                f"        \t shell's relative thickness = {geometry['thickness']}\n"
                f"        \t power = {geometry.get('power', 2.0)}"
            )
        else:
            density_logic = (
                f"         density type = {geometry['type']} ;\n"
                f"        \t Y = {geometry['thickness']}"
            )

        opacity_logic = (f"        \t tau(min) = {dust['opacity']} ; "
                         f"tau(max) = {dust['opacity']}")

        return {
            "{{SPECTRAL_LOGIC}}": spectral_logic,
            "{{BB_LOGIC}}": bb_logic,
            "{{TEMPERATURE_LOGIC}}": temp_logic,
            "{{LUMINOSITY_LOGIC}}": lum_logic,
            "{{ABSORPTION_LOGIC}}": abs_logic,
            "{{DUST_TEMP_LOGIC}}": f"                            Td = {dust['temperature']} K",
            "{{OPTICAL_PROPERTIES_LOGIC}}": (
                f"        Optical properties index = "
                f"{dust.get('properties_index', 'common_and_addl_grain_composite')}"
            ),
            "{{COMPOSITION_LOGIC}}": comp_logic,
            "{{ABUNDANCES_LOGIC}}": abund_logic,
            "{{SIZE_DIST_LOGIC}}": sd_logic,
            "{{SIZE_PARAMS_LOGIC}}": sd_params,
            "{{TSUB_LOGIC}}": f"        Tsub = {dust.get('sublimation_temp', 1500.0)}",
            "{{DENSITY_LOGIC}}": density_logic,
            "{{OPACITY_LOGIC}}": opacity_logic
        }