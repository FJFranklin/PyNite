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
    """
#%%
    __steel = None
    __aluminium = None
#%%
    @staticmethod
    def steel():
        """
        Default material: steel (E=209GPa; G=79.3GPa)

        Returns
        -------
        material : Material
            object with material properties for steel
        """
        if Material.__steel is None:
            Material.__steel = Material(210E9, 79.3E9)
        return Material.__steel

#%%
    @staticmethod
    def aluminium():
        """
        Default material: aluminium (E=70GPa; G=26GPa)

        Returns
        -------
        material : Material
            object with material properties for aluminium
        """
        if Material.__aluminium is None:
            Material.__aluminium = Material(70E9, 26E9)
        return Material.__aluminium

#%%
    def __init__(self, E, G):
        """
        Parameters
        ----------
        E : number
            Young's elastic modulus
        G : number
            shear modulus
        """
        self.E = E
        self.G = G
