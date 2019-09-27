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
        Torsion constant (rather than related polar second moment of area) about x-axis
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
            view.Add(triangle, [color])

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

#%%
    def _torsional_stress(self, Mx, y, z):
        # This only really works for circular sections
        tau_xy = -z * Mx / self.J
        tau_zx =  y * Mx / self.J

        return tau_xy, tau_zx

#%%
    def _shear_stress(self, material, Fy, Fz, y, z):
        tau_xy = 0
        tau_zx = 0
        return tau_xy, tau_zx

#%%
    def _stress(self, XFM, result, material):
        x  = XFM[0]
        Fx = XFM[1]
        Fy = XFM[2]
        Fz = XFM[3]
        Mx = XFM[4]
        My = XFM[5]
        Mz = XFM[6]

        count = len(self._outline[0])

        v = np.zeros((count,3))
        s = np.zeros(count)

        for i in range(0, count):
            coord = self._outline[0][i]
            y = float(coord[0])
            z = float(coord[1])

            sigma_x = Fx / self.A - y * Mz / self.Izz - z * My / self.Iyy

            tau_xy, tau_zx = self._shear_stress(material, Fy, Fz, y, z)

            txy, tzx = self._torsional_stress(Mx, y, z)
            tau_xy += txy
            tau_zx += tzx

            if result == 'seq':
                s[i] = (sigma_x**2 + 3 * (tau_xy**2 + tau_zx**2))**0.5
            elif result == 'sxx':
                s[i] = sigma_x
            elif result == 'txy':
                s[i] = tau_xy
            else: # if result == 'tzx':
                s[i] = tau_zx

            v[i,0] = x
            v[i,1] = y
            v[i,2] = z

        if result == 'seq':
            c_min = 0
            c_max = material.sigma_y
        elif result == 'sxx':
            c_min = -material.sigma_y
            c_max =  material.sigma_y
        else: # use shear yield stress for normalisation
            c_min = -material.sigma_y / 1.73205080757
            c_max =  material.sigma_y / 1.73205080757

        s = (np.minimum(c_max, np.maximum(c_min, s)) - c_min) / (c_max - c_min)

        return v, s

#%%
    def DisplayResults(self, view, basis, forces_moments, result, material):
        """
        Displays the section in 3D indicating stress

        Parameters
        ----------
        view : Viewer3D
            3D viewer for plotting members
        basis : MemberBasis
            Member3D or basis object with coordinate system and member length
        forces_moments : [x,Fx,Fy,Fz,Mx,My,Mz] array
            (Nx,7) array of shear forces and bending moments
        result : str
            Stress to display - one of: 'seq' (von Mises equivalent), 'sxx' (axial)
        material : Material
            Material to use for section
        """

        if self._outline is None: # not enough information
            self.__line(view, basis.ToGlobal(np.asarray([[0,0,0],[basis.L,0,0]])), [0,1])
        else:
            vertex = np.zeros((3,3))
            vcolor = np.zeros((3,4))

            v0, s0 = self._stress(forces_moments[0], result, material)

            v0 = basis.ToGlobal(v0)
            s0 = view.ColorFromValue(s0)

            xcount = len(forces_moments)
            icount = len(v0)

            for x in range(1, xcount):
                vx, sx = self._stress(forces_moments[x], result, material)

                vx = basis.ToGlobal(vx)
                sx = view.ColorFromValue(sx)

                for i in range(0, icount):
                    if i == 0:
                        j = icount - 1
                    else:
                        j = i - 1

                    vertex[0,:] = v0[i,:]
                    vertex[1,:] = v0[j,:]
                    vertex[2,:] = vx[j,:]
                    vcolor[0,:] = s0[i,:]
                    vcolor[1,:] = s0[j,:]
                    vcolor[2,:] = sx[j,:]
                    view.Add(vertex, vcolor)

                    vertex[0,:] = v0[i,:]
                    vertex[1,:] = vx[j,:]
                    vertex[2,:] = vx[i,:]
                    vcolor[0,:] = s0[i,:]
                    vcolor[1,:] = sx[j,:]
                    vcolor[2,:] = sx[i,:]
                    view.Add(vertex, vcolor)

                v0, s0 = vx, sx

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
    def outline(b, d, r, t=0):
        hw = b / 2
        ht = d / 2

        if r > 0:
            rw = hw - r
            rt = ht - r
            sr = r / 2
            cr = r * 0.86602540378

            pts = [ [-rw, ht], [-rw-sr, rt+cr], [-rw-cr, rt+sr], [-hw, rt], [-hw,-rt], [-rw-cr,-rt-sr], [-rw-sr,-rt-cr], [-rw,-ht],
                    [ rw,-ht], [ rw+sr,-rt-cr], [ rw+cr,-rt-sr], [ hw,-rt], [ hw, rt], [ rw+cr, rt+sr], [ rw+sr, rt+cr], [ rw, ht]  ]
        elif t > 0:
            pts = [ [hw,ht], [hw-t,ht], [0,ht], [-hw+t,ht], [-hw,ht], [-hw,ht-t], [-hw,0], [-hw,-ht+t], [-hw,-ht], [-hw+t,-ht], [0,-ht], [hw-t,-ht], [hw,-ht], [hw,-ht+t], [hw,0], [hw,ht-t] ]
        else:
            pts = [ [hw,ht], [0,ht], [-hw,ht], [-hw,0], [-hw,-ht], [0,-ht], [hw,-ht], [hw,0] ]

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

        a = max([breadth, depth]) / 2
        b = min([breadth, depth]) / 2
        e = b / a

        Iyy, Izz, A = Rectangular.properties(breadth, depth, radius)

        self.Iyy = Iyy
        self.Izz = Izz
        self.J   = a * b**3 * (16/3 - 3.36 * e * (1 - e**4 / 12))
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
        self.J   = 2 * (breadth + depth) * thickness**3 / 3
        self.A   = out_A - in_A

        self._outline = [ Rectangular.outline(breadth, depth, radius, thickness), Rectangular.outline(breadth - 2 * thickness, depth - 2 * thickness, radius - thickness) ]

        self.breadth   = breadth
        self.depth     = depth
        self.radius    = radius
        self.thickness = thickness

#%%
    def _shear_stress(self, material, Fy, Fz, y, z):
        tau_xy = 0
        tau_zx = 0

        sign = -1

        if y < 0:
            y = -y
            sign = -sign
        if z < 0:
            z = -z
            sign = -sign

        if self.radius == 0:
            # Stress at surface, estimated from FEA + thin-walled shear-flow theory
            # TODO: check signs

            z_tau_xy = sign * (Fz / self.Iyy) * ((self.depth - self.thickness) / 2) * y
            z_tau_zx = (Fz / self.Iyy) * (((self.depth - self.thickness) / 2) * (self.breadth / 2 - self.thickness) + ((self.depth / 2 + z) / 2) * (self.depth / 2 - z))

            y_tau_xy = (Fy / self.Izz) * (((self.breadth - self.thickness) / 2) * (self.depth / 2 - self.thickness) + ((self.breadth / 2 + y) / 2) * (self.breadth / 2 - y))
            y_tau_zx = sign * (Fy / self.Izz) * ((self.breadth - self.thickness) / 2) * z

            # corner stress effect:
            eb = (self.breadth / 2 - y) / self.thickness
            ed = (self.depth / 2   - z) / self.thickness
            if eb < 2.2:
                ebFz = eb * (1 / 1.1 - eb / 4.84)
            else:
                ebFz = 1
            if eb < 2:
                ebFx = eb * (1 - eb / 4)
            else:
                ebFx = 1
            if ed < 2:
                edFz = ed * (1 - ed / 4)
            else:
                edFz = 1
            if ed < 2.2:
                edFx = ed * (1 / 1.1 - ed / 4.84)
            else:
                edFx = 1

            tau_xy = z_tau_xy * ebFz + y_tau_xy * edFx
            tau_zx = z_tau_zx * edFz + y_tau_zx * ebFx

        return tau_xy, tau_zx

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

#%%
    def _shear_stress(self, material, Fy, Fz, y, z):
        r2 = y * y + z * z
        z2 = z * z
        y2 = y * y
        v = material.v
        # From Timoshenko & Goodier, Theory of Elasticity
        # TODO: check signs
        y_tau_xy = (-Fy / self.Izz) * ((1 + 2 * v) / (4 * (1 + v))) * z * y
        y_tau_zx = ( Fy / self.Izz) * ((3 + 2 * v) / (8 * (1 + v))) * (r2 - y2 - z2 * (1 - 2 * v) / (3 + 2 * v))
        z_tau_xy = (-Fz / self.Iyy) * ((1 + 2 * v) / (4 * (1 + v))) * z * y
        z_tau_zx = ( Fz / self.Iyy) * ((3 + 2 * v) / (8 * (1 + v))) * (r2 - z2 - y2 * (1 - 2 * v) / (3 + 2 * v))
        return (y_tau_xy + z_tau_xy), (y_tau_zx + z_tau_zx)

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
        self.J   = np.pi * diameter * thickness**3 / 3
        self.A   = out_A - in_A

        self._outline = [ Circular.outline(diameter), Circular.outline(diameter - 2 * thickness) ]

        self.diameter  = diameter
        self.thickness = thickness

#%%
    def _shear_stress(self, material, Fy, Fz, y, z):
        r2 = y * y + z * z
        s2 = z * z / r2
        c2 = y * y / r2
        sc = y * z / r2
        # Formula estimated from FEA results; effect of thickness not included
        # TODO: check signs
        y_mean = Fy / self.A
        z_mean = Fz / self.A
        tau_xy = y_mean * s2 * 2 - z_mean * sc * 2
        tau_zx = z_mean * c2 * 2 - y_mean * sc * 2
        return tau_xy, tau_zx

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
        self.J   = (2 * breadth * flange**3 + (depth - flange) * web**3) / 3
        self.A   = out_A - in_A

        self._outline = [ Universal.outline(breadth, depth, flange, web, radius) ]

        self.breadth = breadth
        self.depth   = depth
        self.flange  = flange
        self.web     = web
        self.radius  = radius
