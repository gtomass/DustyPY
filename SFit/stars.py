class Star():

    def __init__(self, Name = '', Temperature = 0, Luminosity = 0):
        self._Name = Name
        self._Temperature = Temperature
        self._Luminosity = Luminosity
    
    def __str__(self):
        return fr'{self._Name}: Temperature = {self._Temperature} K and Luminosity = {self._Luminosity} L_sol'
    
    def set_Name(self, Name):
        self._Name = Name

    def get_Name(self):
        return self._Name

    def set_Temperature(self, Temperature):
        self._Temperature = Temperature
    
    def get_Temperature(self):
        return self._Temperature

    def set_Luminosity(self, Luminosity):
        self._Luminosity = Luminosity

    def get_Luminosity(self):
        return self._Luminosity


class Dust():

    def __init__(self, tau = 0.1, DustSize = {}, Temperature = 800, Composition = {}):
        self._tau = tau
        self._DustSize = DustSize
        self._Temperature = Temperature
        self._Composition = Composition

    def set_tau(self, tau):
        self._tau = tau
    
    def get_tau(self):
        return self._tau

    def set_DustSize(self, DustSize):
        self._DustSize = DustSize
    
    def get_DustSize(self):
        return self._DustSize

    def set_Temperature(self, Temperature):
        self._Temperature = Temperature
    
    def get_Temperature(self):
        return self._Temperature

    def set_Composition(self, Composition):
        self._Composition = Composition

    def get_Composition(self):
        return self._Composition


class Model():

    def __init__(self, Name = '', NbStar = 0, Stars = list(), Dust = Dust(), distance = 1):
        self._Name = Name
        self._NbStar = NbStar
        self._Stars = Stars
        self._Dust = Dust
        self._distance = distance


        self.__Check()

    def set_Name(self, Name):
        self._Name = Name

    def get_Name(self):
        return self._Name
    
    def set_NbStar(self, NBStar):
        self._NbStar = NBStar
        self.__Check()

    def get_NbStar(self):
        return self._NbStar

    def set_Stars(self,Stars = list()):
        self._Stars = Stars
        self.__Check()
    
    def get_Stars(self):
        return self._Stars
    
    def set_Dust(self,Dust = Dust()):
        self._Dust = Dust
    
    def get_Dust(self):
        return self._Dust
    
    def set_Distance(self, distance):
        self._distance = distance
    
    def get_Distance(self):
        return self._distance
    
    def __Check(self):

        if len(self._Stars) != self._NbStar:
            raise ValueError("Number of stars don't match")
        if self._Dust.get_tau() < 0:
            raise ValueError("Tau must be positive")
        if self._distance < 0:
            raise ValueError("Distance must be positive")
        


if __name__== '__main__':
    S = Star()

    S.set_Name('AFGL 4106')
    print(S.get_Name())