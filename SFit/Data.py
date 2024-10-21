from . import utils as utils

class Data():

    def __init__(self, xdata = None, ydata = None, xerr = None, yerr = None):
        self._xdata = xdata
        self._ydata = ydata
        self._xerr = xerr
        self._yerr = yerr

    def set_xdata(self, xdata):
        self._xdata = xdata

    def set_ydata(self, ydata):
        self._ydata = ydata

    def set_xerr(self, xerr):
        self._xerr = xerr

    def set_yerr(self, yerr):
        self._yerr = yerr

    def get_xdata(self):
        return self._xdata
    
    def get_ydata(self):
        return self._ydata
    
    def get_xerr(self):
        return self._xerr
    
    def get_yerr(self):
        return self._yerr
    
    def ImportData(self, Path):
        extention = Path.split('.')[-1]

        if extention == 'csv':
            return utils.LoadCSV(Path)
        elif extention in ['fits', 'fit']:
            return utils.LoadFits(Path)
        else:
            return utils.LoadFile(Path)

    def ScatterData(self, unit=None, xlim=None, ylim=None, ax=None, scale='linear', kwargs=None):
        if self._yerr is None:
            utils.ScatterPlot(self._ydata,self._xdata,unit=unit,xlim=xlim,ylim=ylim,ax=ax,scale=scale,kwargs=kwargs)
        else:
            print('implementer errorbar')
            utils.ScatterPlot(self._ydata,self._xdata,unit=unit,xlim=xlim,ylim=ylim,ax=ax,scale=scale,kwargs=kwargs)

    def UnredData(self, EBV, Rv = 3.1):
        self._ydata = utils.Unred(self._xdata*1e4,self._ydata, EBV=EBV,Rv=Rv)
        self._yerr = utils.Unred(self._xdata*1e4,self._yerr, EBV=EBV,Rv=Rv) if self._yerr is not None else None

    def ConvertJansky(self):
        self._ydata = utils.WattToJansky(self._ydata,self._xdata*1e-6)
        self._yerr = utils.WattToJansky(self._yerr,self._xdata*1e-6) if self._yerr is not None else None
        # try:
        #     self._ydata = utils.WattToJansky(self._ydata.value,self._xdata.value)
        # except Exception:
            
    def ConvertWatt(self):
        self._ydata = utils.JanskyToWatt(self._ydata,self._xdata*1e-6)
        self._yerr = utils.JanskyToWatt(self._yerr,self._xdata*1e-6) if self._yerr is not None else None
        
        # try:
        #     self._ydata = utils.WattToJansky(self._ydata.value,self._xdata.value)
        # except Exception: