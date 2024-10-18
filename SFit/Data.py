class Data():

    def __init__(self, xdata = None, ydata = None):
        self._xdata = xdata
        self._ydata = ydata

    def set_xdata(self, xdata):
        self._xdata = xdata

    def set_ydata(self, ydata):
        self._ydata = ydata

    def get_xdata(self):
        return self._xdata
    
    def get_ydata(self):
        return self._ydata
