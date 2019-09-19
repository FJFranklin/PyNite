# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 2019

@author: Francis James Franklin
"""

import numpy as np

# %%
class Section(object):
    """
    Base class for describing sections

    Attributes
    ----------
    Iyy : number
        Second moment of area for bending about y-axis
    Izz : number
        Second moment of area for bending about z-axis
    J : number
        Polar second moment of area for torsion about x-axis
    A : number
        Cross-sectional area
    """
#%%
    def __init__(self):
        self.Iyy = 0
        self.Izz = 0
        self.J   = 0
        self.A   = 0

        self._outline = None

#%%  
    def __line(self, view, vertices, indices):
        if len(indices) == 2: # just draw a single line between two points
            view.Line(vertices[indices])
        else: # close the loop, and draw
            indices = [*indices, indices[0]]
            view.Line(vertices[indices])

#%%  
    def __face(self, view, vertices, indices, color, center=None):
        # first close the loop
        indices = [*indices, indices[0]]

        if center is None:
            center = vertices[indices[0]]
            vertex = vertices[indices[1:]]
            vcount = len(indices) - 1
        else:
            vertex = vertices[indices]
            vcount = len(indices)

        for i in range(0, vcount-1):
            triangle = np.asarray([center,vertex[i],vertex[i+1]])
            view.Add(triangle, color)

#%%  
    def Display(self, view, wireframe, basis, color):
        """
        Displays the section in 3D

        Parameters
        ----------
        view : Viewer3D
            3D viewer for plotting members
        wireframe : boolean
            If true, only plot wireframe
        basis : MemberBasis
            Member3D or basis object with coordinate system and member length
        color : [r,g,b,a] array
            Color to use for displaying section
        """

        if self._outline is None:
            self.__line(view, basis.ToGlobal(np.asarray([[0,0,0],[basis.L,0,0]])), [0,1])
        else:
            outline = self._outline[0]
            count   = len(outline)

            end_pts = np.zeros((count*2, 3))

            for i in range(0, count):
                p = outline[i]
                end_pts[i,1:3]       = p[:]
                end_pts[count+i,1:3] = p[:]
                end_pts[count+i,0]   = basis.L

            end_pts = basis.ToGlobal(end_pts)

            # Tube

            for i in range(0, count):
                self.__line(view, end_pts, [i,i+count])

                if not wireframe:
                    if i == 0:
                        indices = [count-1, count+count-1, count, 0]
                    else:
                        indices = [i-1, count+i-1, count+i, i]

                    self.__face(view, end_pts, indices, color)

            # Ends

            self.__line(view, end_pts, range(0,count))
            self.__line(view, end_pts, range(2*count-1,count-1,-1))

            if len(self._outline) > 1: # Hollow section
                outline = self._outline[1]
                count   = len(outline)

                end_pts = np.zeros((count*2, 3))

                for i in range(0, count):
                    p = outline[i]
                    end_pts[i,1:3]       = p[:]
                    end_pts[count+i,1:3] = p[:]
                    end_pts[count+i,0]   = basis.L

                end_pts = basis.ToGlobal(end_pts)

                # Ends

                self.__line(view, end_pts, range(0,count))
                self.__line(view, end_pts, range(2*count-1,count-1,-1))

# %%
def Generic(Iyy, Izz, J, A):
    """
    Generic section with only the essential properties

    Returns
    -------
    section : Section
        object with section properties
    """
    S = Section()
    S.Iyy = Iyy
    S.Izz = Izz
    S.J   = J
    S.A   = A
    return S

# %%
class Rectangular(Section):
    """
    Solid rectangular section

    Attributes
    ----------
    breadth : number
        breadth in local y-direction
    depth : number
        depth in local z-direction
    radius : number
        corner radius [default: 0]
    """
#%%
    @staticmethod
    def outline(b, d, r):
        hw = b / 2
        ht = d / 2

        if r > 0:
            rw = hw - r
            rt = ht - r
            sr = r / 2
            cr = r * 0.86602540378

            pts = [ [-rw, ht], [-rw-sr, rt+cr], [-rw-cr, rt+sr], [-hw, rt], [-hw,-rt], [-rw-cr,-rt-sr], [-rw-sr,-rt-cr], [-rw,-ht],
                    [ rw,-ht], [ rw+sr,-rt-cr], [ rw+cr,-rt-sr], [ hw,-rt], [ hw, rt], [ rw+cr, rt+sr], [ rw+sr, rt+cr], [ rw, ht]  ]
        else:
            pts = [ [hw,ht], [-hw,ht], [-hw,-ht], [hw,-ht] ]

        return pts

    @staticmethod
    def properties(b, d, r):
        hw = b / 2
        ht = d / 2

        Iyy = b * d**3 / 12
        Izz = b**3 * d / 12
        A   = b * d

        if r > 0:
            common = r**4 * (np.pi/4 - 16/(9*np.pi) - 1/3)
            Iyy += common - r**2 * ((d - r)**2 - np.pi * (ht - r + 4*r/(3*np.pi))**2)
            Izz += common - r**2 * ((b - r)**2 - np.pi * (hw - r + 4*r/(3*np.pi))**2)
            A   -= (4 - np.pi) * r**2

        return Iyy, Izz, A

#%%
    def __init__(self, breadth, depth, radius=0):
        """
        Parameters
        ----------
        breadth : number
            breadth in local y-direction
        depth : number
            depth in local z-direction
        radius : number
            outer radius of section corners, optional
        """
        Section.__init__(self)

        Iyy, Izz, A = Rectangular.properties(breadth, depth, radius)

        self.Iyy = Iyy
        self.Izz = Izz
        self.J   = self.Iyy + self.Izz
        self.A   = A

        self._outline = [ Rectangular.outline(breadth, depth, radius) ]

        self.breadth = breadth
        self.depth   = depth
        self.radius  = radius

# %%
class RHS(Section):
    """
    Rectangular hollow section

    Attributes
    ----------
    breadth : number
        breadth in local y-direction
    depth : number
        depth in local z-direction
    thickness : number
        thickness of section
    radius : number
        corner radius [default: 0]
    """
#%%
    def __init__(self, breadth, depth, thickness, radius=0):
        """
        Parameters
        ----------
        breadth : number
            breadth in local y-direction
        depth : number
            depth in local z-direction
        thickness : number
            thickness of section
        radius : number
            outer radius of section corners, optional
        """
        Section.__init__(self)

        out_Iyy, out_Izz, out_A = Rectangular.properties(breadth, depth, radius)
        in_Iyy,  in_Izz,  in_A  = Rectangular.properties(breadth - 2 * thickness, depth - 2 * thickness, radius - thickness)

        self.Iyy = out_Iyy - in_Iyy
        self.Izz = out_Izz - in_Izz
        self.J   = self.Iyy + self.Izz
        self.A   = out_A - in_A

        self._outline = [ Rectangular.outline(breadth, depth, radius), Rectangular.outline(breadth - 2 * thickness, depth - 2 * thickness, radius - thickness) ]

        self.breadth   = breadth
        self.depth     = depth
        self.radius    = radius
        self.thickness = thickness

# %%
class Circular(Section):
    """
    Solid circular section

    Attributes
    ----------
    diameter : number
        diameter
    """
#%%
    @staticmethod
    def outline(d):
        hw = d / 2
        ht = d / 2

        rw = 0
        rt = 0
        sr = ht / 2
        cr = ht * 0.86602540378

        pts = [ [-rw, ht], [-rw-sr, rt+cr], [-rw-cr, rt+sr], [-hw, rt], [-rw-cr,-rt-sr], [-rw-sr,-rt-cr],
                [ rw,-ht], [ rw+sr,-rt-cr], [ rw+cr,-rt-sr], [ hw,-rt], [ rw+cr, rt+sr], [ rw+sr, rt+cr]  ]

        return pts

    @staticmethod
    def properties(d):
        Iyy = np.pi * d**4 / 64
        Izz = Iyy
        A   = np.pi * d**2 / 4

        return Iyy, Izz, A

#%%
    def __init__(self, diameter):
        """
        Parameters
        ----------
        diameter : number
            diameter
        """
        Section.__init__(self)

        Iyy, Izz, A = Circular.properties(diameter)

        self.Iyy = Iyy
        self.Izz = Izz
        self.J   = self.Iyy + self.Izz
        self.A   = A

        self._outline = [ Circular.outline(diameter) ]

        self.diameter = diameter

# %%
class CHS(Section):
    """
    Circular hollow section

    Attributes
    ----------
    diameter : number
        outer diameter
    thickness : number
        thickness of section
    """
#%%
    def __init__(self, diameter, thickness):
        """
        Parameters
        ----------
        diameter : number
            outer diameter
        thickness : number
            thickness of section
        """
        Section.__init__(self)

        out_Iyy, out_Izz, out_A = Circular.properties(diameter)
        in_Iyy,  in_Izz,  in_A  = Circular.properties(diameter - 2 * thickness)

        self.Iyy = out_Iyy - in_Iyy
        self.Izz = out_Izz - in_Izz
        self.J   = self.Iyy + self.Izz
        self.A   = out_A - in_A

        self._outline = [ Circular.outline(diameter), Circular.outline(diameter - 2 * thickness) ]

        self.diameter  = diameter
        self.thickness = thickness

# %%
class Universal(Section):
    """
    Rectangular hollow section

    Attributes
    ----------
    breadth : number
        breadth in local y-direction
    depth : number
        depth in local z-direction
    flange : number
        thickness of flange
    web : number
        thickness of web
    radius : number
        corner radius [default: 0]
    """
#%%
    @staticmethod
    def outline(b, d, t, s, r):
        hb = b / 2
        hd = d / 2
        hs = s / 2

        if r > 0:
            ob = hs + r
            od = hd - t - r
            sr = r / 2
            cr = r * 0.86602540378

            pts = [ [-hb, hd], [-hb, hd-t], [-ob, hd-t], [-ob+sr, od+cr], [-ob+cr, od+sr], [-hs, od], [-hs,-od], [-ob+cr,-od-sr], [-ob+sr,-od-cr], [-ob,-hd+t], [-hb,-hd+t], [-hb,-hd],
                    [ hb,-hd], [ hb,-hd+t], [ ob,-hd+t], [ ob-sr,-od-cr], [ ob-cr,-od-sr], [ hs,-od], [ hs, od], [ ob-cr, od+sr], [ ob-sr, od+cr], [ ob, hd-t], [ hb, hd-t], [ hb, hd] ]
        else:
            pts = [ [-hb, hd], [-hb, hd-t], [-hs, hd-t], [-hs,-hd+t], [-hb,-hd+t], [-hb, hd],
                    [ hb,-hd], [ hb,-hd+t], [ hs,-hd+t], [ hs, hd-t], [ hb, hd-t], [ hb,-hd] ]

        return pts

#%%
    def __init__(self, breadth, depth, flange, web, radius=0):
        """
        Parameters
        ----------
        breadth : number
            breadth in local y-direction
        depth : number
            depth in local z-direction
        flange : number
            thickness of flange
        web : number
            thickness of web
        radius : number
            root radius, optional
        """
        Section.__init__(self)

        out_Iyy, out_Izz, out_A = Rectangular.properties(breadth, depth, radius)
        in_Iyy,  in_Izz,  in_A  = Rectangular.properties(breadth - web, depth - 2 * flange, radius)

        self.Iyy = out_Iyy - in_Iyy
        self.Izz = out_Izz - in_Izz
        self.J   = self.Iyy + self.Izz
        self.A   = out_A - in_A

        self._outline = [ Universal.outline(breadth, depth, flange, web, radius) ]

        self.breadth = breadth
        self.depth   = depth
        self.flange  = flange
        self.web     = web
        self.radius  = radius
