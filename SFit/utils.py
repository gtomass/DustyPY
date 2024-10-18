import matplotlib.pyplot as plt
import numpy as np
import subprocess
import glob
import pandas as pd
import astropy
import os



def ScatterPlot(Flux, Wavelength, unit = {'x':'Wavelength', 'y': 'Flux'}, scale='linear',kwargs={}):
    """Create a scatter plot

    Args:
        Flux (array): Flux 
        Wavelength (array): Wavelength
        unit (dict, optional): Description of what is show on each axes. Defaults to {'x':'Wavelength', 'y': 'Flux'}.
        scale (str, optional): Y axe scale. Defaults to 'linear'.
        kwargs (dict, optional): matplolib.pyplot kwargs. Defaults to {}.
    """
    fig,ax = plt.subplots()

    ax.scatter(Wavelength,Flux,**kwargs)
    ax.set_xlabel(f'{unit['x']}')
    ax.set_ylabel(f'{unit['y']}')

    ax.set_yscale(scale)

    plt.show()

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

def LoadFile(file):
    """Load a file in an array

    Args:
        file (string): Path to the file 

    Returns:
        array: array containing the lines of the file
    """

    with open(file , 'r') as f:
        lines = f.readlines()
        return lines
    
def LoadCSV(file, sep=','):
    """Load a csv file

    Args:
        file (string): Path to the file 

    Returns:
        Dataframe: dataframe containing the csv file
    """
    return pd.read_csv(file, sep=sep)

def LoadFits(file):
    """Load a fits file

    Args:
        file (string): Path to the file 

    Returns:
        array: array containing the file
    """
    return astropy.io.fits.open(file)[0].data 

def SearchLine(file,line):
    for i,lines in enumerate(file):
        if line in lines:
            return i 
    raise Exception(f'This line does dot exist: {line}')

# def SearchStr(line,car):
#     for i,seq in enumerate(line.split(' ')):
#         if car in seq:
#             return i
#     return None

# def ModifyStr(line, index, car):
#     l = line.split(' ')
#     l[index] = car
#     return ' '.join(l)

def SaveFile(Path,file):
    with open(Path , 'w') as f:
            f.write("".join(file))

def ChangeParameter(Path,change,car):
     
    file = LoadFile(Path)
    for param in change.keys():
        line = SearchLine(file,car[param])
        new_line = change[param]
        file[line] = new_line

    SaveFile(Path,file)

if __name__=="__main__":
     pass
    #  f = LoadFile('/Users/gabriel/Documents/SFit/AFGL4106.inp')
    #  i = SearchLine(f, 'Number of additional components')
    #  j = SearchStr(f[i], '.nk')  
    #  print(f[i])
    #  print(ModifyStr(f[i], j, 'prout'))