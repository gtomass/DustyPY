class Star():
    """
    Class representing a star with a name, temperature, and luminosity.
    """

    def __init__(self, Name='', Temperature=0, Luminosity=0):
        """
        Initializes an instance of the Star class.

        Parameters:
        Name (str, optional): The name of the star. Defaults to an empty string.
        Temperature (float, optional): The temperature of the star in Kelvin. Defaults to 0.
        Luminosity (float, optional): The luminosity of the star in solar units. Defaults to 0.
        """
        self._Name = Name
        self._Temperature = Temperature
        self._Luminosity = Luminosity
    
    def __str__(self):
        """
        Returns a string representation of the star.

        Returns:
        str: A string describing the star with its name, temperature, and luminosity.
        """
        return fr'{self._Name}: Temperature = {self._Temperature} K and Luminosity = {self._Luminosity} L_sol'
    
    def set_Name(self, Name):
        """
        Sets the name of the star.

        Parameters:
        Name (str): The name of the star.
        """
        self._Name = Name

    def get_Name(self):
        """
        Returns the name of the star.

        Returns:
        str: The name of the star.
        """
        return self._Name

    def set_Temperature(self, Temperature):
        """
        Sets the temperature of the star.

        Parameters:
        Temperature (float): The temperature of the star in Kelvin.
        """
        self._Temperature = Temperature
    
    def get_Temperature(self):
        """
        Returns the temperature of the star.

        Returns:
        float: The temperature of the star in Kelvin.
        """
        return self._Temperature

    def set_Luminosity(self, Luminosity):
        """
        Sets the luminosity of the star.

        Parameters:
        Luminosity (float): The luminosity of the star in solar units.
        """
        self._Luminosity = Luminosity

    def get_Luminosity(self):
        """
        Returns the luminosity of the star.

        Returns:
        float: The luminosity of the star in solar units.
        """
        return self._Luminosity


class Dust():
    """
    Class representing dust with opacity, size, temperature, and composition.
    """

    def __init__(self, tau=0.1, DustSize={}, Temperature=800, Composition={}):
        """
        Initializes an instance of the Dust class.

        Parameters:
        tau (float, optional): The opacity of the dust. Defaults to 0.1.
        DustSize (dict, optional): A dictionary representing the size of the dust. Defaults to an empty dictionary.
        Temperature (float, optional): The temperature of the dust in Kelvin. Defaults to 800.
        Composition (dict, optional): A dictionary representing the composition of the dust. Defaults to an empty dictionary.
        """
        self._tau = tau
        self._DustSize = DustSize
        self._Temperature = Temperature
        self._Composition = Composition

    def set_tau(self, tau):
        """
        Sets the opacity (tau) of the dust.

        Parameters:
        tau (float): The opacity of the dust.
        """
        self._tau = tau

    def get_tau(self):
        """
        Returns the opacity (tau) of the dust.

        Returns:
        float: The opacity of the dust.
        """
        return self._tau

    def set_DustSize(self, DustSize):
        """
        Sets the size of the dust.

        Parameters:
        DustSize (dict): A dictionary representing the size of the dust.
        """
        self._DustSize = DustSize

    def get_DustSize(self):
        """
        Returns the size of the dust.

        Returns:
        dict: A dictionary representing the size of the dust.
        """
        return self._DustSize

    def set_Temperature(self, Temperature):
        """
        Sets the temperature of the dust.

        Parameters:
        Temperature (float): The temperature of the dust in Kelvin.
        """
        self._Temperature = Temperature

    def get_Temperature(self):
        """
        Returns the temperature of the dust.

        Returns:
        float: The temperature of the dust in Kelvin.
        """
        return self._Temperature

    def set_Composition(self, Composition):
        """
        Sets the composition of the dust.

        Parameters:
        Composition (dict): A dictionary representing the composition of the dust.
        """
        self._Composition = Composition

    def get_Composition(self):
        """
        Returns the composition of the dust.

        Returns:
        dict: A dictionary representing the composition of the dust.
        """
        return self._Composition


class Model():
    """
    Class representing a model with a name, number of stars, stars, dust, and distance.
    """

    def __init__(self, Name='', NbStar=0, Stars=list(), Dust=Dust(), distance=1):
        """
        Initializes an instance of the Model class.

        Parameters:
        Name (str, optional): The name of the model. Defaults to an empty string.
        NbStar (int, optional): The number of stars in the model. Defaults to 0.
        Stars (list, optional): A list of stars in the model. Defaults to an empty list.
        Dust (Dust, optional): An instance of the Dust class representing the dust in the model. Defaults to a new Dust instance.
        distance (float, optional): The distance of the model in arbitrary units. Defaults to 1.
        """
        self._Name = Name
        self._NbStar = NbStar
        self._Stars = Stars
        self._Dust = Dust
        self._distance = distance
        self.__Check()

    def set_Name(self, Name):
        """
        Sets the name of the model.

        Parameters:
        Name (str): The name of the model.
        """
        self._Name = Name

    def get_Name(self):
        """
        Returns the name of the model.

        Returns:
        str: The name of the model.
        """
        return self._Name
    
    def set_NbStar(self, NbStar):
        """
        Sets the number of stars in the model.

        Parameters:
        NbStar (int): The number of stars in the model.
        """
        self._NbStar = NbStar
        self.__Check()

    def get_NbStar(self):
        """
        Returns the number of stars in the model.

        Returns:
        int: The number of stars in the model.
        """
        return self._NbStar

    def set_Stars(self, Stars=list()):
        """
        Sets the list of stars in the model.

        Parameters:
        Stars (list, optional): A list of stars in the model. Defaults to an empty list.
        """
        self._Stars = Stars
        self.__Check()
    
    def get_Stars(self):
        """
        Returns the list of stars in the model.

        Returns:
        list: The list of stars in the model.
        """
        return self._Stars
    
    def set_Dust(self, Dust=Dust()):
        """
        Sets the dust in the model.

        Parameters:
        Dust (Dust, optional): An instance of the Dust class representing the dust in the model. Defaults to a new Dust instance.
        """
        self._Dust = Dust
    
    def get_Dust(self):
        """
        Returns the dust in the model.

        Returns:
        Dust: An instance of the Dust class representing the dust in the model.
        """
        return self._Dust
    
    def set_Distance(self, distance):
        """
        Sets the distance of the model.

        Parameters:
        distance (float): The distance of the model in arbitrary units.
        """
        self._distance = distance
    
    def get_Distance(self):
        """
        Returns the distance of the model.

        Returns:
        float: The distance of the model in arbitrary units.
        """
        return self._distance
    
    def __Check(self):
        """
        Checks the consistency of the model's attributes.
        """
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