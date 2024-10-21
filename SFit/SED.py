import SFit.utils as utils

class SED():

    def __init__(self, Flux = list, wavelength = list):
        self._Flux = Flux
        self._Wavelength = wavelength
        
    def set_Flux(self, Flux):
        self._Flux = Flux

    def set_Wavelength(self, wavelength):
        self._Wavelength = wavelength

    def get_Flux(self):
        return self._Flux
    
    def get_Wavelength(self):
        return self._Wavelength

    def PlotSED(self, unit=None, xlim=None, ylim=None, ax=None, scale='linear', kwargs=None):
        utils.Plot(self._Flux,self._Wavelength,unit=unit,xlim=xlim,ylim=ylim,ax=ax,scale=scale,kwargs=kwargs)

if __name__=='__main__':
    f= [1,2,3,4]
    w = [0.1,0.2,0.4,0.8]

    k = {'marker':'+', 'color':'red'}

    s = SED(Flux=f,wavelength=w)
    s.PlotSED(kwargs=k)