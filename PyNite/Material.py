# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 2019

@author: Francis James Franklin
"""
# %%
class Material(object):
    """
    Class object for holding essential material properties.

    Attributes
    ----------
    E : number
        Young's elastic modulus
    G : number
        shear modulus
    v : number
        Poisson's ratio
    density : number
        density
    sigma_y : number
        yield stress [default: 200E6]
    """
#%%
    __steel = None
    __aluminium = None
#%%
    @staticmethod
    def steel():
        """
        Default material: steel (E=209GPa; G=79.3GPa; density=7850kg/m3; yield=235MPa)

        Returns
        -------
        material : Material
            object with material properties for steel
        """
        if Material.__steel is None:
            Material.__steel = Material(210E9, 79.3E9, 7850)
            Material.__steel.sigma_y = 235E6
        return Material.__steel

#%%
    @staticmethod
    def aluminium():
        """
        Default material: aluminium (E=70GPa; G=26GPa; density=2700kg/m3; yield=320MPa)

        Returns
        -------
        material : Material
            object with material properties for aluminium
        """
        if Material.__aluminium is None:
            Material.__aluminium = Material(70E9, 26E9, 2700)
            Material.__aluminium.set_color((0.9, 0.9, 0.9))
            Material.__aluminium.sigma_y = 320E6
        return Material.__aluminium

#%%
    def __init__(self, E, G, density=0):
        """
        Parameters
        ----------
        E : number
            Young's elastic modulus
        G : number
            shear modulus
        density : number
            density, optional [default is 0]
        """
        self.E = E
        self.G = G

        self.v = E / (2 * G) - 1

        self.density = density
        self.sigma_y = 200E6
        self.__color = [0.6,0.6,0.6,1]

#%%
    def set_color(self, rgb):
        """
        Set the color to be used for displaying the material

        Parameters
        ----------
        rgb : (r,g,b) triple
            color to use for material; 0 <= r,g,b <= 1
        """
        r, g, b = rgb

        self.__color = [r,g,b,1]

#%%
    def color(self):
        """
        Get the color to be used for displaying the material

        Returns
        -------
        rgba : [r,g,b,a] array
            color to use for material; 0 <= r,g,b <= 1, a = 1
        """
        return self.__color
