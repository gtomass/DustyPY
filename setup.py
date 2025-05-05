from setuptools import setup, find_packages, Extension

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
        "DustyPY": ["filter/*", "libs/*"],
    },
    ext_modules=[simpson_extension],  # Include the C extension
)