"""
Setup script for DustyPY.
Handles the compilation of the C Simpson library with OpenMP support.
"""

import os
import sys
import tempfile
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

class BuildExt(build_ext):
    """
    Custom build class to detect and enable OpenMP support across 
    different platforms (Windows, Linux, macOS).
    """
    def build_extensions(self):
        compiler_flags = []
        
        # Define flags based on compiler type and OS
        if self.compiler.compiler_type == 'msvc':
            # Windows / MSVC
            compiler_flags.append((['/openmp'], []))
        else:
            # Linux & others / GCC or Clang
            compiler_flags.append((['-fopenmp'], ['-fopenmp']))
            if sys.platform == 'darwin':
                # Special handling for macOS (Homebrew libomp is usually required)
                compiler_flags.append((['-Xpreprocessor', '-fopenmp'], ['-lomp']))

        supported_flags = None
        for c_flags, l_flags in compiler_flags:
            if self.check_openmp_support(c_flags, l_flags):
                supported_flags = (c_flags, l_flags)
                break
        
        if supported_flags:
            print(f"DEBUG: OpenMP found. Compiling with {supported_flags}")
            for ext in self.extensions:
                ext.extra_compile_args.extend(supported_flags[0])
                ext.extra_link_args.extend(supported_flags[1])
        else:
            print("WARNING: OpenMP not found. Compiling without OpenMP support.")
            
        super().build_extensions()

    def check_openmp_support(self, compile_args, link_args):
        """Checks if the compiler can successfully build a basic OpenMP program."""
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                test_file = os.path.join(tmpdir, 'test_openmp.c')
                with open(test_file, 'w', encoding='utf-8') as f:
                    f.write(
                        "#include <omp.h>\n"
                        "#include <stdio.h>\n"
                        "int main() { return 0; }\n"
                    )
                objects = self.compiler.compile(
                    [test_file], 
                    output_dir=tmpdir, 
                    extra_postargs=compile_args
                )
                self.compiler.link_executable(
                    objects, 
                    os.path.join(tmpdir, 'test_openmp'), 
                    extra_postargs=link_args
                )
            return True
        except Exception:
            return False

# Define the C extension for Simpson's rule
simpson_extension = Extension(
    name="dustypy.libs.simpson",
    sources=["src/dustypy/libs/simpson.c"], # Matches your package structure
    include_dirs=["src/dustypy/libs"]
)

setup(
    # Metadata and dependencies are handled by pyproject.toml
    # but we define the C extension and custom build command here
    ext_modules=[simpson_extension],
    cmdclass={'build_ext': BuildExt},
    packages=find_packages(),
    zip_safe=False,
)