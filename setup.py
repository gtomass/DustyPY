import os
import sys
import tempfile
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

class BuildExt(build_ext):
    def build_extensions(self):
        compiler_flags = []
        if self.compiler.compiler_type == 'msvc':
            compiler_flags.append((['/openmp'], []))
        else:
            compiler_flags.append((['-fopenmp'], ['-fopenmp']))
            if sys.platform == 'darwin':
                compiler_flags.append((['-Xpreprocessor', '-fopenmp'], ['-lomp']))

        supported_flags = None
        for c_flags, l_flags in compiler_flags:
            if self.check_openmp_support(c_flags, l_flags):
                supported_flags = (c_flags, l_flags)
                break
        
        if supported_flags:
            print(f"OpenMP found. Compiling with {supported_flags}")
            for ext in self.extensions:
                ext.extra_compile_args.extend(supported_flags[0])
                ext.extra_link_args.extend(supported_flags[1])
        else:
            print("WARNING: OpenMP not found. Compiling without OpenMP support.")
        super().build_extensions()

    def check_openmp_support(self, compile_args, link_args):
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                test_file = os.path.join(tmpdir, 'test_openmp.c')
                with open(test_file, 'w') as f:
                    f.write(
                        "#include <omp.h>\n"
                        "#include <stdio.h>\n"
                        "int main() { return 0; }\n"
                    )
                objects = self.compiler.compile([test_file], output_dir=tmpdir, extra_postargs=compile_args)
                self.compiler.link_executable(objects, os.path.join(tmpdir, 'test_openmp'), extra_postargs=link_args)
            return True
        except Exception:
            return False

# Define the C extension
simpson_extension = Extension(
    name="DustyPY.libs.simpson",  # Module name
    sources=["DustyPY/libs/simpson.c"],  # Path to the C source file
)

# Setup configuration
setup(
    name="DustyPY",
    version="0.1.0",
    author="Gabriel TOMASSINI",
    author_email="gabriel.tomassini@oca.eu",
    description="DustyPY is a convenient Package to run and fit SED using the radiative transfer modeling code dusty",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/gtomass/DustyPY.git",
    package_dir={"DustyPY": "DustyPY"},
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.12',
    install_requires=[
        "matplotlib",
        "numpy",
        "pandas",
        "astropy",
        "scipy",
        "pymcmcstat",
        "synphot"
    ],
    include_package_data=True,
    package_data={
        "DustyPY": ["filter/*", "libs/simpson.c"],
    },
    ext_modules=[simpson_extension],  # Include the C extension
    cmdclass={'build_ext': BuildExt},
)