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
    def __init__(self, **kwargs):
        self.Iyy = 0
        self.Izz = 0
        self.J   = 0
        self.A   = 0

        if 'Iyy' in kwargs:
            self.Iyy = kwargs['Iyy']
        if 'Izz' in kwargs:
            self.Izz = kwargs['Izz']
        if 'A' in kwargs:
            self.A = kwargs['A']
        if 'J' in kwargs:
            self.J = kwargs['J']

        self._outline = None

#%%
    def _set_or_compare(self, Iyy, Izz, A, J, print_comparison=False):
        if self.Iyy > 0:
            err = 100 * (Iyy - self.Iyy) / self.Iyy
            if print_comparison and abs(err) > 0.5:
                print('Iyy: spec={s:.3e}, calc={c:.3e}, diff={d:.1f}%'.format(s=self.Iyy, c=Iyy, d=err))
        else:
            self.Iyy = Iyy

        if self.Izz > 0:
            err = 100 * (Izz - self.Izz) / self.Izz
            if print_comparison and abs(err) > 0.5:
                print('Izz: spec={s:.3e}, calc={c:.3e}, diff={d:.1f}%'.format(s=self.Izz, c=Izz, d=err))
        else:
            self.Izz = Izz

        if self.A > 0:
            err = 100 * (A - self.A) / self.A
            if print_comparison and abs(err) > 0.5:
                print('A:   spec={s:.3e}, calc={c:.3e}, diff={d:.1f}%'.format(s=self.A, c=A, d=err))
        else:
            self.A = A

        if self.J > 0:
            err = 100 * (J - self.J) / self.J
            if print_comparison and abs(err) > 0.5:
                print('J:   spec={s:.3e}, calc={c:.3e}, diff={d:.1f}%'.format(s=self.J, c=J, d=err))
        else:
            self.J = J

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

            #print('{x},{y},{z},{fy},{fz},{sxx},{txy},{tzx}'.format(x=x, y=y, z=z, fy=Fy, fz=Fz, sxx=sigma_x, txy=tau_xy, tzx=tau_zx))

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
    def __init__(self, breadth, depth, radius=0, **kwargs):
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
        Section.__init__(self, **kwargs)

        a = max([breadth, depth]) / 2
        b = min([breadth, depth]) / 2
        e = b / a

        Iyy, Izz, A = Rectangular.properties(breadth, depth, radius)

        J = a * b**3 * (16/3 - 3.36 * e * (1 - e**4 / 12))

        self._set_or_compare(Iyy, Izz, A, J)

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
    def __init__(self, breadth, depth, thickness, radius=0, **kwargs):
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
        Section.__init__(self, **kwargs)

        out_Iyy, out_Izz, out_A = Rectangular.properties(breadth, depth, radius)
        in_Iyy,  in_Izz,  in_A  = Rectangular.properties(breadth - 2 * thickness, depth - 2 * thickness, radius - thickness)

        perimeter = 2 * (breadth + depth)
        if radius > thickness:
            perimeter += 2 * np.pi * (radius - thickness / 2) - 8 * radius

        self._set_or_compare(out_Iyy - in_Iyy, out_Izz - in_Izz, out_A - in_A, (out_A + in_A)**2 * thickness / perimeter)

        self._outline = [ Rectangular.outline(breadth, depth, radius, thickness), Rectangular.outline(breadth - 2 * thickness, depth - 2 * thickness, radius - thickness) ]

        self.breadth   = breadth
        self.depth     = depth
        self.radius    = radius
        self.thickness = thickness

    def __shear(self, Iyy, b, d, Fz, y, z, sign): # shear stresses (at surface only)
        t = self.thickness
        r = self.radius

        if r > t: # rounded (hollow) rectangle - guesswork... TODO: Check/Correct!
            ylim = (b / 2 - r)
            zlim = (d / 2 - r)
            if y < ylim: # top edge
                tau_xy = sign * (Fz / Iyy) * ((d - t) / 2) * y
                tau_zx = 0
            elif z < zlim: # right edge
                a1t_top    = ((d - t) / 2) * ylim
                a1t_corner = (np.pi / 4) * zlim * (2 * r - t) + (r * (r - t) + t * t / 3)
                a1t_right  = ((zlim + z) / 2) * (zlim - z)
                tau_xy = 0
                tau_zx = (Fz / Iyy) * (a1t_top + a1t_corner + a1t_right)
            else:
                a1t_top    = ((d - t) / 2) * ylim
                a1t_corner = (np.pi / 4) * zlim * (2 * r - t) + (r * (r - t) + t * t / 3)
                lim_xy = sign * (Fz / Iyy) * a1t_top
                lim_zx = (Fz / Iyy) * (a1t_top + a1t_corner)
                tau_xy = lim_xy * (z - zlim) / r
                tau_zx = lim_zx * (y - ylim) / r

        else: # non-rounded RHS
            # Stress at surface, estimated from FEA + thin-walled shear-flow theory
            # TODO: check signs
            tau_xy = sign * (Fz / Iyy) * ((d - t) / 2) * y
            tau_zx = (Fz / Iyy) * (((d - t) / 2) * (b / 2 - t) + ((d / 2 + z) / 2) * (d / 2 - z))

            # corner stress effect:
            eb = (b / 2 - y) / t
            ed = (d / 2 - z) / t

            if eb < 2.2:
                tau_xy *= eb * (1 / 1.1 - eb / 4.84)
            if ed < 2:
                tau_zx *= ed * (1 - ed / 4)

        return tau_xy, tau_zx

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

        z_tau_xy, z_tau_zx = self.__shear(self.Iyy, self.breadth, self.depth, Fz, y, z, sign)
        y_tau_zx, y_tau_xy = self.__shear(self.Izz, self.depth, self.breadth, Fy, z, y, sign)

        return (z_tau_xy + y_tau_xy), (z_tau_zx + y_tau_zx)

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
    def __init__(self, diameter, **kwargs):
        """
        Parameters
        ----------
        diameter : number
            diameter
        """
        Section.__init__(self, **kwargs)

        Iyy, Izz, A = Circular.properties(diameter)

        self._set_or_compare(Iyy, Izz, A, Iyy + Izz)

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
    def __init__(self, diameter, thickness, **kwargs):
        """
        Parameters
        ----------
        diameter : number
            outer diameter
        thickness : number
            thickness of section
        """
        Section.__init__(self, **kwargs)

        out_Iyy, out_Izz, out_A = Circular.properties(diameter)
        in_Iyy,  in_Izz,  in_A  = Circular.properties(diameter - 2 * thickness)

        self._set_or_compare(out_Iyy - in_Iyy, out_Izz - in_Izz, out_A - in_A, out_Iyy - in_Iyy + out_Izz - in_Izz)

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
        y_mean = -Fy / self.A
        z_mean = -Fz / self.A
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
    def __init__(self, breadth, depth, flange, web, radius=0, **kwargs):
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
        Section.__init__(self, **kwargs)

        out_Iyy, out_Izz, out_A = Rectangular.properties(breadth, depth, 0)
        in_Iyy,  in_,     in_A  = Rectangular.properties(breadth - web, depth - 2 * flange, radius)

        in_1,    in_Izz,  in_2  = Rectangular.properties(breadth, depth - 2 * flange, 0)
        in_1,   web_Izz,  in_2  = Rectangular.properties(web, depth - 2 * flange, 0)

        self._set_or_compare(out_Iyy - in_Iyy, out_Izz - in_Izz + web_Izz, out_A - in_A, (2 * breadth * flange**3 + (depth - flange) * web**3) / 3)

        self._outline = [ Universal.outline(breadth, depth, flange, web, radius) ]

        self.breadth = breadth
        self.depth   = depth
        self.flange  = flange
        self.web     = web
        self.radius  = radius

class BlueBook(object):
    '''A selection of EuroCode 3 sections'''

    __UniversalBeam = None

    @staticmethod
    def __create_UniversalBeam():
        S = {}
        S['1016 x 305'] = {}
        S['1016 x 305']['x 584'] = Universal(3.140e-01, 1.056e+00, 6.400e-02, 3.600e-02, 3.000e-02, area=7.440e-02, Iyy=1.246e-02, Izz=3.340e-04, J=7.150e-05)
        S['1016 x 305']['x 494'] = Universal(3.090e-01, 1.036e+00, 5.400e-02, 3.100e-02, 3.000e-02, area=6.290e-02, Iyy=1.028e-02, Izz=2.680e-04, J=4.400e-05)
        S['1016 x 305']['x 438'] = Universal(3.050e-01, 1.026e+00, 4.900e-02, 2.690e-02, 3.000e-02, area=5.560e-02, Iyy=9.100e-03, Izz=2.340e-04, J=3.190e-05)
        S['1016 x 305']['x 415'] = Universal(3.040e-01, 1.020e+00, 4.600e-02, 2.600e-02, 3.000e-02, area=5.290e-02, Iyy=8.530e-03, Izz=2.170e-04, J=2.700e-05)
        S['1016 x 305']['x 393'] = Universal(3.030e-01, 1.016e+00, 4.390e-02, 2.440e-02, 3.000e-02, area=5.000e-02, Iyy=8.080e-03, Izz=2.050e-04, J=2.330e-05)
        S['1016 x 305']['x 350'] = Universal(3.020e-01, 1.008e+00, 4.000e-02, 2.110e-02, 3.000e-02, area=4.450e-02, Iyy=7.230e-03, Izz=1.850e-04, J=1.720e-05)
        S['1016 x 305']['x 314'] = Universal(3.000e-01, 9.999e-01, 3.590e-02, 1.910e-02, 3.000e-02, area=4.000e-02, Iyy=6.440e-03, Izz=1.620e-04, J=1.260e-05)
        S['1016 x 305']['x 272'] = Universal(3.000e-01, 9.901e-01, 3.100e-02, 1.650e-02, 3.000e-02, area=3.470e-02, Iyy=5.540e-03, Izz=1.400e-04, J=8.350e-06)
        S['1016 x 305']['x 249'] = Universal(3.000e-01, 9.801e-01, 2.600e-02, 1.650e-02, 3.000e-02, area=3.170e-02, Iyy=4.810e-03, Izz=1.180e-04, J=5.820e-06)
        S['1016 x 305']['x 222'] = Universal(3.000e-01, 9.703e-01, 2.110e-02, 1.600e-02, 3.000e-02, area=2.830e-02, Iyy=4.080e-03, Izz=9.550e-05, J=3.900e-06)
        S['914 x 419'] = {}
        S['914 x 419']['x 388'] = Universal(4.205e-01, 9.210e-01, 3.660e-02, 2.140e-02, 2.410e-02, area=4.940e-02, Iyy=7.200e-03, Izz=4.540e-04, J=1.730e-05)
        S['914 x 419']['x 343'] = Universal(4.185e-01, 9.118e-01, 3.200e-02, 1.940e-02, 2.410e-02, area=4.370e-02, Iyy=6.260e-03, Izz=3.920e-04, J=1.190e-05)
        S['914 x 305'] = {}
        S['914 x 305']['x 576'] = Universal(3.220e-01, 9.930e-01, 6.500e-02, 3.610e-02, 1.900e-02, area=7.330e-02, Iyy=1.102e-02, Izz=3.650e-04, J=7.130e-05)
        S['914 x 305']['x 521'] = Universal(3.190e-01, 9.810e-01, 5.890e-02, 3.300e-02, 1.900e-02, area=6.640e-02, Iyy=9.820e-03, Izz=3.210e-04, J=5.340e-05)
        S['914 x 305']['x 474'] = Universal(3.160e-01, 9.710e-01, 5.410e-02, 3.000e-02, 1.900e-02, area=6.040e-02, Iyy=8.860e-03, Izz=2.870e-04, J=4.100e-05)
        S['914 x 305']['x 425'] = Universal(3.130e-01, 9.610e-01, 4.900e-02, 2.690e-02, 1.900e-02, area=5.420e-02, Iyy=7.880e-03, Izz=2.520e-04, J=3.030e-05)
        S['914 x 305']['x 381'] = Universal(3.100e-01, 9.510e-01, 4.390e-02, 2.440e-02, 1.900e-02, area=4.860e-02, Iyy=6.970e-03, Izz=2.190e-04, J=2.200e-05)
        S['914 x 305']['x 345'] = Universal(3.080e-01, 9.430e-01, 3.990e-02, 2.210e-02, 1.900e-02, area=4.400e-02, Iyy=6.260e-03, Izz=1.950e-04, J=1.650e-05)
        S['914 x 305']['x 313'] = Universal(3.090e-01, 9.320e-01, 3.450e-02, 2.110e-02, 1.900e-02, area=3.980e-02, Iyy=5.480e-03, Izz=1.700e-04, J=1.160e-05)
        S['914 x 305']['x 289'] = Universal(3.077e-01, 9.266e-01, 3.200e-02, 1.950e-02, 1.910e-02, area=3.680e-02, Iyy=5.040e-03, Izz=1.560e-04, J=9.260e-06)
        S['914 x 305']['x 271'] = Universal(3.070e-01, 9.230e-01, 3.000e-02, 1.840e-02, 1.900e-02, area=3.460e-02, Iyy=4.720e-03, Izz=1.450e-04, J=7.690e-06)
        S['914 x 305']['x 253'] = Universal(3.055e-01, 9.184e-01, 2.790e-02, 1.730e-02, 1.910e-02, area=3.230e-02, Iyy=4.360e-03, Izz=1.330e-04, J=6.260e-06)
        S['914 x 305']['x 238'] = Universal(3.050e-01, 9.150e-01, 2.590e-02, 1.650e-02, 1.900e-02, area=3.040e-02, Iyy=4.060e-03, Izz=1.230e-04, J=5.140e-06)
        S['914 x 305']['x 224'] = Universal(3.041e-01, 9.104e-01, 2.390e-02, 1.590e-02, 1.910e-02, area=2.860e-02, Iyy=3.760e-03, Izz=1.120e-04, J=4.220e-06)
        S['914 x 305']['x 201'] = Universal(3.033e-01, 9.030e-01, 2.020e-02, 1.510e-02, 1.910e-02, area=2.560e-02, Iyy=3.250e-03, Izz=9.420e-05, J=2.910e-06)
        S['838 x 292'] = {}
        S['838 x 292']['x 226'] = Universal(2.938e-01, 8.509e-01, 2.680e-02, 1.610e-02, 1.780e-02, area=2.890e-02, Iyy=3.400e-03, Izz=1.140e-04, J=5.140e-06)
        S['838 x 292']['x 194'] = Universal(2.924e-01, 8.407e-01, 2.170e-02, 1.470e-02, 1.780e-02, area=2.470e-02, Iyy=2.790e-03, Izz=9.070e-05, J=3.060e-06)
        S['838 x 292']['x 176'] = Universal(2.917e-01, 8.349e-01, 1.880e-02, 1.400e-02, 1.780e-02, area=2.240e-02, Iyy=2.460e-03, Izz=7.800e-05, J=2.210e-06)
        S['762 x 267'] = {}
        S['762 x 267']['x 197'] = Universal(2.680e-01, 7.698e-01, 2.540e-02, 1.560e-02, 1.650e-02, area=2.510e-02, Iyy=2.400e-03, Izz=8.170e-05, J=4.040e-06)
        S['762 x 267']['x 173'] = Universal(2.667e-01, 7.622e-01, 2.160e-02, 1.430e-02, 1.650e-02, area=2.200e-02, Iyy=2.050e-03, Izz=6.850e-05, J=2.670e-06)
        S['762 x 267']['x 147'] = Universal(2.652e-01, 7.540e-01, 1.750e-02, 1.280e-02, 1.650e-02, area=1.870e-02, Iyy=1.690e-03, Izz=5.460e-05, J=1.590e-06)
        S['762 x 267']['x 134'] = Universal(2.644e-01, 7.500e-01, 1.550e-02, 1.200e-02, 1.650e-02, area=1.710e-02, Iyy=1.510e-03, Izz=4.790e-05, J=1.190e-06)
        S['686 x 254'] = {}
        S['686 x 254']['x 170'] = Universal(2.558e-01, 6.929e-01, 2.370e-02, 1.450e-02, 1.520e-02, area=2.170e-02, Iyy=1.700e-03, Izz=6.630e-05, J=3.080e-06)
        S['686 x 254']['x 152'] = Universal(2.545e-01, 6.875e-01, 2.100e-02, 1.320e-02, 1.520e-02, area=1.940e-02, Iyy=1.500e-03, Izz=5.780e-05, J=2.200e-06)
        S['686 x 254']['x 140'] = Universal(2.537e-01, 6.835e-01, 1.900e-02, 1.240e-02, 1.520e-02, area=1.780e-02, Iyy=1.360e-03, Izz=5.180e-05, J=1.690e-06)
        S['686 x 254']['x 125'] = Universal(2.530e-01, 6.779e-01, 1.620e-02, 1.170e-02, 1.520e-02, area=1.590e-02, Iyy=1.180e-03, Izz=4.380e-05, J=1.160e-06)
        S['610 x 305'] = {}
        S['610 x 305']['x 238'] = Universal(3.114e-01, 6.358e-01, 3.140e-02, 1.840e-02, 1.650e-02, area=3.030e-02, Iyy=2.090e-03, Izz=1.580e-04, J=7.850e-06)
        S['610 x 305']['x 179'] = Universal(3.071e-01, 6.202e-01, 2.360e-02, 1.410e-02, 1.650e-02, area=2.280e-02, Iyy=1.530e-03, Izz=1.140e-04, J=3.400e-06)
        S['610 x 305']['x 149'] = Universal(3.048e-01, 6.124e-01, 1.970e-02, 1.180e-02, 1.650e-02, area=1.900e-02, Iyy=1.260e-03, Izz=9.310e-05, J=2.000e-06)
        S['610 x 229'] = {}
        S['610 x 229']['x 140'] = Universal(2.302e-01, 6.172e-01, 2.210e-02, 1.310e-02, 1.270e-02, area=1.780e-02, Iyy=1.120e-03, Izz=4.510e-05, J=2.160e-06)
        S['610 x 229']['x 125'] = Universal(2.290e-01, 6.122e-01, 1.960e-02, 1.190e-02, 1.270e-02, area=1.590e-02, Iyy=9.860e-04, Izz=3.930e-05, J=1.540e-06)
        S['610 x 229']['x 113'] = Universal(2.282e-01, 6.076e-01, 1.730e-02, 1.110e-02, 1.270e-02, area=1.440e-02, Iyy=8.730e-04, Izz=3.430e-05, J=1.110e-06)
        S['610 x 229']['x 101'] = Universal(2.276e-01, 6.026e-01, 1.480e-02, 1.050e-02, 1.270e-02, area=1.290e-02, Iyy=7.580e-04, Izz=2.910e-05, J=7.700e-07)
        S['610 x 178'] = {}
        S['610 x 178']['x 100'] = Universal(1.792e-01, 6.074e-01, 1.720e-02, 1.130e-02, 1.270e-02, area=1.280e-02, Iyy=7.250e-04, Izz=1.660e-05, J=9.500e-07)
        S['610 x 178']['x 92'] = Universal(1.788e-01, 6.030e-01, 1.500e-02, 1.090e-02, 1.270e-02, area=1.170e-02, Iyy=6.460e-04, Izz=1.440e-05, J=7.100e-07)
        S['610 x 178']['x 82'] = Universal(1.779e-01, 5.986e-01, 1.280e-02, 1.000e-02, 1.270e-02, area=1.040e-02, Iyy=5.590e-04, Izz=1.210e-05, J=4.880e-07)
        S['533 x 312'] = {}
        S['533 x 312']['x 273'] = Universal(3.202e-01, 5.771e-01, 3.760e-02, 2.110e-02, 1.270e-02, area=3.480e-02, Iyy=1.990e-03, Izz=2.060e-04, J=1.290e-05)
        S['533 x 312']['x 219'] = Universal(3.174e-01, 5.603e-01, 2.920e-02, 1.830e-02, 1.270e-02, area=2.790e-02, Iyy=1.510e-03, Izz=1.560e-04, J=6.420e-06)
        S['533 x 312']['x 182'] = Universal(3.145e-01, 5.507e-01, 2.440e-02, 1.520e-02, 1.270e-02, area=2.310e-02, Iyy=1.230e-03, Izz=1.270e-04, J=3.730e-06)
        S['533 x 312']['x 151'] = Universal(3.120e-01, 5.425e-01, 2.030e-02, 1.270e-02, 1.270e-02, area=1.920e-02, Iyy=1.010e-03, Izz=1.030e-04, J=2.160e-06)
        S['533 x 210'] = {}
        S['533 x 210']['x 138'] = Universal(2.139e-01, 5.491e-01, 2.360e-02, 1.470e-02, 1.270e-02, area=1.760e-02, Iyy=8.610e-04, Izz=3.860e-05, J=2.500e-06)
        S['533 x 210']['x 122'] = Universal(2.119e-01, 5.445e-01, 2.130e-02, 1.270e-02, 1.270e-02, area=1.550e-02, Iyy=7.600e-04, Izz=3.390e-05, J=1.780e-06)
        S['533 x 210']['x 109'] = Universal(2.108e-01, 5.395e-01, 1.880e-02, 1.160e-02, 1.270e-02, area=1.390e-02, Iyy=6.680e-04, Izz=2.940e-05, J=1.260e-06)
        S['533 x 210']['x 101'] = Universal(2.100e-01, 5.367e-01, 1.740e-02, 1.080e-02, 1.270e-02, area=1.290e-02, Iyy=6.150e-04, Izz=2.690e-05, J=1.010e-06)
        S['533 x 210']['x 92'] = Universal(2.093e-01, 5.331e-01, 1.560e-02, 1.010e-02, 1.270e-02, area=1.170e-02, Iyy=5.520e-04, Izz=2.390e-05, J=7.570e-07)
        S['533 x 210']['x 82'] = Universal(2.088e-01, 5.283e-01, 1.320e-02, 9.600e-03, 1.270e-02, area=1.050e-02, Iyy=4.750e-04, Izz=2.010e-05, J=5.150e-07)
        S['533 x 165'] = {}
        S['533 x 165']['x 85'] = Universal(1.665e-01, 5.349e-01, 1.650e-02, 1.030e-02, 1.270e-02, area=1.080e-02, Iyy=4.850e-04, Izz=1.270e-05, J=7.380e-07)
        S['533 x 165']['x 75'] = Universal(1.659e-01, 5.291e-01, 1.360e-02, 9.700e-03, 1.270e-02, area=9.520e-03, Iyy=4.110e-04, Izz=1.040e-05, J=4.790e-07)
        S['533 x 165']['x 66'] = Universal(1.651e-01, 5.247e-01, 1.140e-02, 8.900e-03, 1.270e-02, area=8.370e-03, Iyy=3.500e-04, Izz=8.590e-06, J=3.200e-07)
        S['457 x 191'] = {}
        S['457 x 191']['x 161'] = Universal(1.994e-01, 4.920e-01, 3.200e-02, 1.800e-02, 1.020e-02, area=2.060e-02, Iyy=7.980e-04, Izz=4.250e-05, J=5.150e-06)
        S['457 x 191']['x 133'] = Universal(1.967e-01, 4.806e-01, 2.630e-02, 1.530e-02, 1.020e-02, area=1.700e-02, Iyy=6.380e-04, Izz=3.350e-05, J=2.920e-06)
        S['457 x 191']['x 106'] = Universal(1.940e-01, 4.692e-01, 2.060e-02, 1.260e-02, 1.020e-02, area=1.350e-02, Iyy=4.890e-04, Izz=2.510e-05, J=1.460e-06)
        S['457 x 191']['x 98'] = Universal(1.928e-01, 4.672e-01, 1.960e-02, 1.140e-02, 1.020e-02, area=1.250e-02, Iyy=4.570e-04, Izz=2.350e-05, J=1.210e-06)
        S['457 x 191']['x 89'] = Universal(1.919e-01, 4.634e-01, 1.770e-02, 1.050e-02, 1.020e-02, area=1.140e-02, Iyy=4.100e-04, Izz=2.090e-05, J=9.070e-07)
        S['457 x 191']['x 82'] = Universal(1.913e-01, 4.600e-01, 1.600e-02, 9.900e-03, 1.020e-02, area=1.040e-02, Iyy=3.710e-04, Izz=1.870e-05, J=6.920e-07)
        S['457 x 191']['x 74'] = Universal(1.904e-01, 4.570e-01, 1.450e-02, 9.000e-03, 1.020e-02, area=9.460e-03, Iyy=3.330e-04, Izz=1.670e-05, J=5.180e-07)
        S['457 x 191']['x 67'] = Universal(1.899e-01, 4.534e-01, 1.270e-02, 8.500e-03, 1.020e-02, area=8.550e-03, Iyy=2.940e-04, Izz=1.450e-05, J=3.710e-07)
        S['457 x 152'] = {}
        S['457 x 152']['x 82'] = Universal(1.553e-01, 4.658e-01, 1.890e-02, 1.050e-02, 1.020e-02, area=1.050e-02, Iyy=3.660e-04, Izz=1.180e-05, J=8.920e-07)
        S['457 x 152']['x 74'] = Universal(1.544e-01, 4.620e-01, 1.700e-02, 9.600e-03, 1.020e-02, area=9.450e-03, Iyy=3.270e-04, Izz=1.050e-05, J=6.590e-07)
        S['457 x 152']['x 67'] = Universal(1.538e-01, 4.580e-01, 1.500e-02, 9.000e-03, 1.020e-02, area=8.560e-03, Iyy=2.890e-04, Izz=9.130e-06, J=4.770e-07)
        S['457 x 152']['x 60'] = Universal(1.529e-01, 4.546e-01, 1.330e-02, 8.100e-03, 1.020e-02, area=7.620e-03, Iyy=2.550e-04, Izz=7.950e-06, J=3.380e-07)
        S['457 x 152']['x 52'] = Universal(1.524e-01, 4.498e-01, 1.090e-02, 7.600e-03, 1.020e-02, area=6.660e-03, Iyy=2.140e-04, Izz=6.450e-06, J=2.140e-07)
        S['406 x 178'] = {}
        S['406 x 178']['x 85'] = Universal(1.819e-01, 4.172e-01, 1.820e-02, 1.090e-02, 1.020e-02, area=1.090e-02, Iyy=3.170e-04, Izz=1.830e-05, J=9.300e-07)
        S['406 x 178']['x 74'] = Universal(1.795e-01, 4.128e-01, 1.600e-02, 9.500e-03, 1.020e-02, area=9.450e-03, Iyy=2.730e-04, Izz=1.550e-05, J=6.280e-07)
        S['406 x 178']['x 67'] = Universal(1.788e-01, 4.094e-01, 1.430e-02, 8.800e-03, 1.020e-02, area=8.550e-03, Iyy=2.430e-04, Izz=1.360e-05, J=4.610e-07)
        S['406 x 178']['x 60'] = Universal(1.779e-01, 4.064e-01, 1.280e-02, 7.900e-03, 1.020e-02, area=7.650e-03, Iyy=2.160e-04, Izz=1.200e-05, J=3.330e-07)
        S['406 x 178']['x 54'] = Universal(1.777e-01, 4.026e-01, 1.090e-02, 7.700e-03, 1.020e-02, area=6.900e-03, Iyy=1.870e-04, Izz=1.020e-05, J=2.310e-07)
        S['406 x 140'] = {}
        S['406 x 140']['x 53'] = Universal(1.433e-01, 4.066e-01, 1.290e-02, 7.900e-03, 1.020e-02, area=6.790e-03, Iyy=1.830e-04, Izz=6.350e-06, J=2.900e-07)
        S['406 x 140']['x 46'] = Universal(1.422e-01, 4.032e-01, 1.120e-02, 6.800e-03, 1.020e-02, area=5.860e-03, Iyy=1.570e-04, Izz=5.380e-06, J=1.900e-07)
        S['406 x 140']['x 39'] = Universal(1.418e-01, 3.980e-01, 8.600e-03, 6.400e-03, 1.020e-02, area=4.970e-03, Iyy=1.250e-04, Izz=4.100e-06, J=1.070e-07)
        S['356 x 171'] = {}
        S['356 x 171']['x 67'] = Universal(1.732e-01, 3.634e-01, 1.570e-02, 9.100e-03, 1.020e-02, area=8.550e-03, Iyy=1.950e-04, Izz=1.360e-05, J=5.570e-07)
        S['356 x 171']['x 57'] = Universal(1.722e-01, 3.580e-01, 1.300e-02, 8.100e-03, 1.020e-02, area=7.260e-03, Iyy=1.600e-04, Izz=1.110e-05, J=3.340e-07)
        S['356 x 171']['x 51'] = Universal(1.715e-01, 3.550e-01, 1.150e-02, 7.400e-03, 1.020e-02, area=6.490e-03, Iyy=1.410e-04, Izz=9.680e-06, J=2.380e-07)
        S['356 x 171']['x 45'] = Universal(1.711e-01, 3.514e-01, 9.700e-03, 7.000e-03, 1.020e-02, area=5.730e-03, Iyy=1.210e-04, Izz=8.110e-06, J=1.580e-07)
        S['356 x 127'] = {}
        S['356 x 127']['x 39'] = Universal(1.260e-01, 3.534e-01, 1.070e-02, 6.600e-03, 1.020e-02, area=4.980e-03, Iyy=1.020e-04, Izz=3.580e-06, J=1.510e-07)
        S['356 x 127']['x 33'] = Universal(1.254e-01, 3.490e-01, 8.500e-03, 6.000e-03, 1.020e-02, area=4.210e-03, Iyy=8.250e-05, Izz=2.800e-06, J=8.790e-08)
        S['305 x 165'] = {}
        S['305 x 165']['x 54'] = Universal(1.669e-01, 3.104e-01, 1.370e-02, 7.900e-03, 8.900e-03, area=6.880e-03, Iyy=1.170e-04, Izz=1.060e-05, J=3.480e-07)
        S['305 x 165']['x 46'] = Universal(1.657e-01, 3.066e-01, 1.180e-02, 6.700e-03, 8.900e-03, area=5.870e-03, Iyy=9.900e-05, Izz=8.960e-06, J=2.220e-07)
        S['305 x 165']['x 40'] = Universal(1.650e-01, 3.034e-01, 1.020e-02, 6.000e-03, 8.900e-03, area=5.130e-03, Iyy=8.500e-05, Izz=7.640e-06, J=1.470e-07)
        S['305 x 127'] = {}
        S['305 x 127']['x 48'] = Universal(1.253e-01, 3.110e-01, 1.400e-02, 9.000e-03, 8.900e-03, area=6.120e-03, Iyy=9.570e-05, Izz=4.610e-06, J=3.180e-07)
        S['305 x 127']['x 42'] = Universal(1.243e-01, 3.072e-01, 1.210e-02, 8.000e-03, 8.900e-03, area=5.340e-03, Iyy=8.200e-05, Izz=3.890e-06, J=2.110e-07)
        S['305 x 127']['x 37'] = Universal(1.234e-01, 3.044e-01, 1.070e-02, 7.100e-03, 8.900e-03, area=4.720e-03, Iyy=7.170e-05, Izz=3.360e-06, J=1.480e-07)
        S['305 x 102'] = {}
        S['305 x 102']['x 33'] = Universal(1.024e-01, 3.127e-01, 1.080e-02, 6.600e-03, 7.600e-03, area=4.180e-03, Iyy=6.500e-05, Izz=1.940e-06, J=1.220e-07)
        S['305 x 102']['x 28'] = Universal(1.018e-01, 3.087e-01, 8.800e-03, 6.000e-03, 7.600e-03, area=3.590e-03, Iyy=5.370e-05, Izz=1.550e-06, J=7.400e-08)
        S['305 x 102']['x 25'] = Universal(1.016e-01, 3.051e-01, 7.000e-03, 5.800e-03, 7.600e-03, area=3.160e-03, Iyy=4.460e-05, Izz=1.230e-06, J=4.770e-08)
        S['254 x 146'] = {}
        S['254 x 146']['x 43'] = Universal(1.473e-01, 2.596e-01, 1.270e-02, 7.200e-03, 7.600e-03, area=5.480e-03, Iyy=6.540e-05, Izz=6.770e-06, J=2.390e-07)
        S['254 x 146']['x 37'] = Universal(1.464e-01, 2.560e-01, 1.090e-02, 6.300e-03, 7.600e-03, area=4.720e-03, Iyy=5.540e-05, Izz=5.710e-06, J=1.530e-07)
        S['254 x 146']['x 31'] = Universal(1.461e-01, 2.514e-01, 8.600e-03, 6.000e-03, 7.600e-03, area=3.970e-03, Iyy=4.410e-05, Izz=4.480e-06, J=8.550e-08)
        S['254 x 102'] = {}
        S['254 x 102']['x 28'] = Universal(1.022e-01, 2.604e-01, 1.000e-02, 6.300e-03, 7.600e-03, area=3.610e-03, Iyy=4.000e-05, Izz=1.790e-06, J=9.570e-08)
        S['254 x 102']['x 25'] = Universal(1.019e-01, 2.572e-01, 8.400e-03, 6.000e-03, 7.600e-03, area=3.200e-03, Iyy=3.410e-05, Izz=1.490e-06, J=6.420e-08)
        S['254 x 102']['x 22'] = Universal(1.016e-01, 2.540e-01, 6.800e-03, 5.700e-03, 7.600e-03, area=2.800e-03, Iyy=2.840e-05, Izz=1.190e-06, J=4.150e-08)
        S['203 x 133'] = {}
        S['203 x 133']['x 30'] = Universal(1.339e-01, 2.068e-01, 9.600e-03, 6.400e-03, 7.600e-03, area=3.820e-03, Iyy=2.900e-05, Izz=3.850e-06, J=1.030e-07)
        S['203 x 133']['x 25'] = Universal(1.332e-01, 2.032e-01, 7.800e-03, 5.700e-03, 7.600e-03, area=3.200e-03, Iyy=2.340e-05, Izz=3.080e-06, J=5.960e-08)
        S['203 x 102'] = {}
        S['203 x 102']['x 23'] = Universal(1.018e-01, 2.032e-01, 9.300e-03, 5.400e-03, 7.600e-03, area=2.940e-03, Iyy=2.100e-05, Izz=1.640e-06, J=7.020e-08)
        S['178 x 102'] = {}
        S['178 x 102']['x 19'] = Universal(1.012e-01, 1.778e-01, 7.900e-03, 4.800e-03, 7.600e-03, area=2.430e-03, Iyy=1.360e-05, Izz=1.370e-06, J=4.410e-08)
        S['152 x 89'] = {}
        S['152 x 89']['x 16'] = Universal(8.870e-02, 1.524e-01, 7.700e-03, 4.500e-03, 7.600e-03, area=2.030e-03, Iyy=8.340e-06, Izz=8.980e-07, J=3.560e-08)
        S['127 x 76'] = {}
        S['127 x 76']['x 13'] = Universal(7.600e-02, 1.270e-01, 7.600e-03, 4.000e-03, 7.600e-03, area=1.650e-03, Iyy=4.730e-06, Izz=5.570e-07, J=2.850e-08)

        return S # No. data rows = 107

    @staticmethod
    def UniversalBeam():
        if BlueBook.__UniversalBeam is None:
            BlueBook.__UniversalBeam = BlueBook.__create_UniversalBeam()
        return BlueBook.__UniversalBeam

    __CHS_ColdFormed = None

    @staticmethod
    def __create_CHS_ColdFormed():
        S = {}
        S['33.7'] = {}
        S['33.7']['3.0'] = CHS(3.370e-02, 3.000e-03, area=2.890e-04, Iyy=3.440e-08, Izz=3.440e-08, J=6.880e-08)
        S['42.4'] = {}
        S['42.4']['3.0'] = CHS(4.240e-02, 3.000e-03, area=3.710e-04, Iyy=7.250e-08, Izz=7.250e-08, J=1.450e-07)
        S['48.3'] = {}
        S['48.3']['3.0'] = CHS(4.830e-02, 3.000e-03, area=4.270e-04, Iyy=1.100e-07, Izz=1.100e-07, J=2.200e-07)
        S['48.3']['4.0'] = CHS(4.830e-02, 4.000e-03, area=5.570e-04, Iyy=1.380e-07, Izz=1.380e-07, J=2.750e-07)
        S['60.3'] = {}
        S['60.3']['3.0'] = CHS(6.030e-02, 3.000e-03, area=5.400e-04, Iyy=2.220e-07, Izz=2.220e-07, J=4.440e-07)
        S['60.3']['4.0'] = CHS(6.030e-02, 4.000e-03, area=7.070e-04, Iyy=2.820e-07, Izz=2.820e-07, J=5.630e-07)
        S['76.1'] = {}
        S['76.1']['3.0'] = CHS(7.610e-02, 3.000e-03, area=6.890e-04, Iyy=4.610e-07, Izz=4.610e-07, J=9.220e-07)
        S['76.1']['4.0'] = CHS(7.610e-02, 4.000e-03, area=9.060e-04, Iyy=5.910e-07, Izz=5.910e-07, J=1.180e-06)
        S['88.9'] = {}
        S['88.9']['3.0'] = CHS(8.890e-02, 3.000e-03, area=8.100e-04, Iyy=7.480e-07, Izz=7.480e-07, J=1.500e-06)
        S['88.9']['3.5'] = CHS(8.890e-02, 3.500e-03, area=9.390e-04, Iyy=8.570e-07, Izz=8.570e-07, J=1.710e-06)
        S['88.9']['4.0'] = CHS(8.890e-02, 4.000e-03, area=1.070e-03, Iyy=9.630e-07, Izz=9.630e-07, J=1.930e-06)
        S['88.9']['5.0'] = CHS(8.890e-02, 5.000e-03, area=1.320e-03, Iyy=1.160e-06, Izz=1.160e-06, J=2.330e-06)
        S['88.9']['6.3'] = CHS(8.890e-02, 6.300e-03, area=1.630e-03, Iyy=1.400e-06, Izz=1.400e-06, J=2.800e-06)
        S['114.3'] = {}
        S['114.3']['3.0'] = CHS(1.143e-01, 3.000e-03, area=1.050e-03, Iyy=1.630e-06, Izz=1.630e-06, J=3.250e-06)
        S['114.3']['3.5'] = CHS(1.143e-01, 3.500e-03, area=1.220e-03, Iyy=1.870e-06, Izz=1.870e-06, J=3.740e-06)
        S['114.3']['4.0'] = CHS(1.143e-01, 4.000e-03, area=1.390e-03, Iyy=2.110e-06, Izz=2.110e-06, J=4.220e-06)
        S['114.3']['5.0'] = CHS(1.143e-01, 5.000e-03, area=1.720e-03, Iyy=2.570e-06, Izz=2.570e-06, J=5.140e-06)
        S['114.3']['6.0'] = CHS(1.143e-01, 6.000e-03, area=2.040e-03, Iyy=3.000e-06, Izz=3.000e-06, J=6.000e-06)
        S['114.3']['6.3'] = CHS(1.143e-01, 6.300e-03, area=2.140e-03, Iyy=3.130e-06, Izz=3.130e-06, J=6.250e-06)
        S['139.7'] = {}
        S['139.7']['3.0'] = CHS(1.397e-01, 3.000e-03, area=1.290e-03, Iyy=3.010e-06, Izz=3.010e-06, J=6.020e-06)
        S['139.7']['4.0'] = CHS(1.397e-01, 4.000e-03, area=1.710e-03, Iyy=3.930e-06, Izz=3.930e-06, J=7.860e-06)
        S['139.7']['5.0'] = CHS(1.397e-01, 5.000e-03, area=2.120e-03, Iyy=4.810e-06, Izz=4.810e-06, J=9.610e-06)
        S['139.7']['6.0'] = CHS(1.397e-01, 6.000e-03, area=2.520e-03, Iyy=5.640e-06, Izz=5.640e-06, J=1.130e-05)
        S['139.7']['6.3'] = CHS(1.397e-01, 6.300e-03, area=2.640e-03, Iyy=5.890e-06, Izz=5.890e-06, J=1.180e-05)
        S['139.7']['8.0'] = CHS(1.397e-01, 8.000e-03, area=3.310e-03, Iyy=7.200e-06, Izz=7.200e-06, J=1.440e-05)
        S['139.7']['10.0'] = CHS(1.397e-01, 1.000e-02, area=4.070e-03, Iyy=8.620e-06, Izz=8.620e-06, J=1.720e-05)
        S['168.3'] = {}
        S['168.3']['4.0'] = CHS(1.683e-01, 4.000e-03, area=2.060e-03, Iyy=6.970e-06, Izz=6.970e-06, J=1.390e-05)
        S['168.3']['4.5'] = CHS(1.683e-01, 4.500e-03, area=2.320e-03, Iyy=7.770e-06, Izz=7.770e-06, J=1.550e-05)
        S['168.3']['5.0'] = CHS(1.683e-01, 5.000e-03, area=2.570e-03, Iyy=8.560e-06, Izz=8.560e-06, J=1.710e-05)
        S['168.3']['6.0'] = CHS(1.683e-01, 6.000e-03, area=3.060e-03, Iyy=1.010e-05, Izz=1.010e-05, J=2.020e-05)
        S['168.3']['6.3'] = CHS(1.683e-01, 6.300e-03, area=3.210e-03, Iyy=1.050e-05, Izz=1.050e-05, J=2.110e-05)
        S['168.3']['8.0'] = CHS(1.683e-01, 8.000e-03, area=4.030e-03, Iyy=1.300e-05, Izz=1.300e-05, J=2.600e-05)
        S['168.3']['10.0'] = CHS(1.683e-01, 1.000e-02, area=4.970e-03, Iyy=1.560e-05, Izz=1.560e-05, J=3.130e-05)
        S['168.3']['12.5'] = CHS(1.683e-01, 1.250e-02, area=6.120e-03, Iyy=1.870e-05, Izz=1.870e-05, J=3.740e-05)
        S['193.7'] = {}
        S['193.7']['4.0'] = CHS(1.937e-01, 4.000e-03, area=2.380e-03, Iyy=1.070e-05, Izz=1.070e-05, J=2.150e-05)
        S['193.7']['4.5'] = CHS(1.937e-01, 4.500e-03, area=2.670e-03, Iyy=1.200e-05, Izz=1.200e-05, J=2.400e-05)
        S['193.7']['5.0'] = CHS(1.937e-01, 5.000e-03, area=2.960e-03, Iyy=1.320e-05, Izz=1.320e-05, J=2.640e-05)
        S['193.7']['6.0'] = CHS(1.937e-01, 6.000e-03, area=3.540e-03, Iyy=1.560e-05, Izz=1.560e-05, J=3.120e-05)
        S['193.7']['6.3'] = CHS(1.937e-01, 6.300e-03, area=3.710e-03, Iyy=1.630e-05, Izz=1.630e-05, J=3.260e-05)
        S['193.7']['8.0'] = CHS(1.937e-01, 8.000e-03, area=4.670e-03, Iyy=2.020e-05, Izz=2.020e-05, J=4.030e-05)
        S['193.7']['10.0'] = CHS(1.937e-01, 1.000e-02, area=5.770e-03, Iyy=2.440e-05, Izz=2.440e-05, J=4.880e-05)
        S['193.7']['12.5'] = CHS(1.937e-01, 1.250e-02, area=7.120e-03, Iyy=2.930e-05, Izz=2.930e-05, J=5.870e-05)
        S['219.1'] = {}
        S['219.1']['4.5'] = CHS(2.191e-01, 4.500e-03, area=3.030e-03, Iyy=1.750e-05, Izz=1.750e-05, J=3.490e-05)
        S['219.1']['5.0'] = CHS(2.191e-01, 5.000e-03, area=3.360e-03, Iyy=1.930e-05, Izz=1.930e-05, J=3.860e-05)
        S['219.1']['6.0'] = CHS(2.191e-01, 6.000e-03, area=4.020e-03, Iyy=2.280e-05, Izz=2.280e-05, J=4.560e-05)
        S['219.1']['6.3'] = CHS(2.191e-01, 6.300e-03, area=4.210e-03, Iyy=2.390e-05, Izz=2.390e-05, J=4.770e-05)
        S['219.1']['8.0'] = CHS(2.191e-01, 8.000e-03, area=5.310e-03, Iyy=2.960e-05, Izz=2.960e-05, J=5.920e-05)
        S['219.1']['10.0'] = CHS(2.191e-01, 1.000e-02, area=6.570e-03, Iyy=3.600e-05, Izz=3.600e-05, J=7.200e-05)
        S['219.1']['12.0'] = CHS(2.191e-01, 1.200e-02, area=7.810e-03, Iyy=4.200e-05, Izz=4.200e-05, J=8.400e-05)
        S['219.1']['12.5'] = CHS(2.191e-01, 1.250e-02, area=8.110e-03, Iyy=4.340e-05, Izz=4.340e-05, J=8.690e-05)
        S['219.1']['16.0'] = CHS(2.191e-01, 1.600e-02, area=1.020e-02, Iyy=5.300e-05, Izz=5.300e-05, J=1.060e-04)
        S['244.5'] = {}
        S['244.5']['5.0'] = CHS(2.445e-01, 5.000e-03, area=3.760e-03, Iyy=2.700e-05, Izz=2.700e-05, J=5.400e-05)
        S['244.5']['6.0'] = CHS(2.445e-01, 6.000e-03, area=4.500e-03, Iyy=3.200e-05, Izz=3.200e-05, J=6.400e-05)
        S['244.5']['6.3'] = CHS(2.445e-01, 6.300e-03, area=4.710e-03, Iyy=3.350e-05, Izz=3.350e-05, J=6.690e-05)
        S['244.5']['8.0'] = CHS(2.445e-01, 8.000e-03, area=5.940e-03, Iyy=4.160e-05, Izz=4.160e-05, J=8.320e-05)
        S['244.5']['10.0'] = CHS(2.445e-01, 1.000e-02, area=7.370e-03, Iyy=5.070e-05, Izz=5.070e-05, J=1.010e-04)
        S['244.5']['12.0'] = CHS(2.445e-01, 1.200e-02, area=8.770e-03, Iyy=5.940e-05, Izz=5.940e-05, J=1.190e-04)
        S['244.5']['12.5'] = CHS(2.445e-01, 1.250e-02, area=9.110e-03, Iyy=6.150e-05, Izz=6.150e-05, J=1.230e-04)
        S['244.5']['16.0'] = CHS(2.445e-01, 1.600e-02, area=1.150e-02, Iyy=7.530e-05, Izz=7.530e-05, J=1.510e-04)
        S['273.0'] = {}
        S['273.0']['4.0'] = CHS(2.730e-01, 4.000e-03, area=3.380e-03, Iyy=3.060e-05, Izz=3.060e-05, J=6.120e-05)
        S['273.0']['4.5'] = CHS(2.730e-01, 4.500e-03, area=3.800e-03, Iyy=3.420e-05, Izz=3.420e-05, J=6.840e-05)
        S['273.0']['5.0'] = CHS(2.730e-01, 5.000e-03, area=4.210e-03, Iyy=3.780e-05, Izz=3.780e-05, J=7.560e-05)
        S['273.0']['6.0'] = CHS(2.730e-01, 6.000e-03, area=5.030e-03, Iyy=4.490e-05, Izz=4.490e-05, J=8.970e-05)
        S['273.0']['6.3'] = CHS(2.730e-01, 6.300e-03, area=5.280e-03, Iyy=4.700e-05, Izz=4.700e-05, J=9.390e-05)
        S['273.0']['8.0'] = CHS(2.730e-01, 8.000e-03, area=6.660e-03, Iyy=5.850e-05, Izz=5.850e-05, J=1.170e-04)
        S['273.0']['10.0'] = CHS(2.730e-01, 1.000e-02, area=8.260e-03, Iyy=7.150e-05, Izz=7.150e-05, J=1.430e-04)
        S['273.0']['12.0'] = CHS(2.730e-01, 1.200e-02, area=9.840e-03, Iyy=8.400e-05, Izz=8.400e-05, J=1.680e-04)
        S['273.0']['12.5'] = CHS(2.730e-01, 1.250e-02, area=1.020e-02, Iyy=8.700e-05, Izz=8.700e-05, J=1.740e-04)
        S['273.0']['16.0'] = CHS(2.730e-01, 1.600e-02, area=1.290e-02, Iyy=1.070e-04, Izz=1.070e-04, J=2.140e-04)
        S['323.9'] = {}
        S['323.9']['5.0'] = CHS(3.239e-01, 5.000e-03, area=5.010e-03, Iyy=6.370e-05, Izz=6.370e-05, J=1.270e-04)
        S['323.9']['6.0'] = CHS(3.239e-01, 6.000e-03, area=5.990e-03, Iyy=7.570e-05, Izz=7.570e-05, J=1.510e-04)
        S['323.9']['6.3'] = CHS(3.239e-01, 6.300e-03, area=6.290e-03, Iyy=7.930e-05, Izz=7.930e-05, J=1.590e-04)
        S['323.9']['8.0'] = CHS(3.239e-01, 8.000e-03, area=7.940e-03, Iyy=9.910e-05, Izz=9.910e-05, J=1.980e-04)
        S['323.9']['10.0'] = CHS(3.239e-01, 1.000e-02, area=9.860e-03, Iyy=1.220e-04, Izz=1.220e-04, J=2.430e-04)
        S['323.9']['12.0'] = CHS(3.239e-01, 1.200e-02, area=1.180e-02, Iyy=1.430e-04, Izz=1.430e-04, J=2.860e-04)
        S['323.9']['12.5'] = CHS(3.239e-01, 1.250e-02, area=1.220e-02, Iyy=1.480e-04, Izz=1.480e-04, J=2.970e-04)
        S['323.9']['16.0'] = CHS(3.239e-01, 1.600e-02, area=1.550e-02, Iyy=1.840e-04, Izz=1.840e-04, J=3.680e-04)
        S['355.6'] = {}
        S['355.6']['5.0'] = CHS(3.556e-01, 5.000e-03, area=5.510e-03, Iyy=8.460e-05, Izz=8.460e-05, J=1.690e-04)
        S['355.6']['6.0'] = CHS(3.556e-01, 6.000e-03, area=6.590e-03, Iyy=1.010e-04, Izz=1.010e-04, J=2.010e-04)
        S['355.6']['6.3'] = CHS(3.556e-01, 6.300e-03, area=6.910e-03, Iyy=1.050e-04, Izz=1.050e-04, J=2.110e-04)
        S['355.6']['8.0'] = CHS(3.556e-01, 8.000e-03, area=8.740e-03, Iyy=1.320e-04, Izz=1.320e-04, J=2.640e-04)
        S['355.6']['10.0'] = CHS(3.556e-01, 1.000e-02, area=1.090e-02, Iyy=1.620e-04, Izz=1.620e-04, J=3.240e-04)
        S['355.6']['12.0'] = CHS(3.556e-01, 1.200e-02, area=1.300e-02, Iyy=1.910e-04, Izz=1.910e-04, J=3.830e-04)
        S['355.6']['12.5'] = CHS(3.556e-01, 1.250e-02, area=1.350e-02, Iyy=1.990e-04, Izz=1.990e-04, J=3.970e-04)
        S['355.6']['16.0'] = CHS(3.556e-01, 1.600e-02, area=1.710e-02, Iyy=2.470e-04, Izz=2.470e-04, J=4.930e-04)
        S['406.4'] = {}
        S['406.4']['6.0'] = CHS(4.064e-01, 6.000e-03, area=7.550e-03, Iyy=1.510e-04, Izz=1.510e-04, J=3.030e-04)
        S['406.4']['6.3'] = CHS(4.064e-01, 6.300e-03, area=7.920e-03, Iyy=1.580e-04, Izz=1.580e-04, J=3.170e-04)
        S['406.4']['8.0'] = CHS(4.064e-01, 8.000e-03, area=1.000e-02, Iyy=1.990e-04, Izz=1.990e-04, J=3.970e-04)
        S['406.4']['10.0'] = CHS(4.064e-01, 1.000e-02, area=1.250e-02, Iyy=2.450e-04, Izz=2.450e-04, J=4.900e-04)
        S['406.4']['12.0'] = CHS(4.064e-01, 1.200e-02, area=1.490e-02, Iyy=2.890e-04, Izz=2.890e-04, J=5.790e-04)
        S['406.4']['12.5'] = CHS(4.064e-01, 1.250e-02, area=1.550e-02, Iyy=3.000e-04, Izz=3.000e-04, J=6.010e-04)
        S['406.4']['16.0'] = CHS(4.064e-01, 1.600e-02, area=1.960e-02, Iyy=3.740e-04, Izz=3.740e-04, J=7.490e-04)
        S['457.0'] = {}
        S['457.0']['6.0'] = CHS(4.570e-01, 6.000e-03, area=8.500e-03, Iyy=2.160e-04, Izz=2.160e-04, J=4.320e-04)
        S['457.0']['6.3'] = CHS(4.570e-01, 6.300e-03, area=8.920e-03, Iyy=2.270e-04, Izz=2.270e-04, J=4.530e-04)
        S['457.0']['8.0'] = CHS(4.570e-01, 8.000e-03, area=1.130e-02, Iyy=2.840e-04, Izz=2.840e-04, J=5.690e-04)
        S['457.0']['10.0'] = CHS(4.570e-01, 1.000e-02, area=1.400e-02, Iyy=3.510e-04, Izz=3.510e-04, J=7.020e-04)
        S['457.0']['12.0'] = CHS(4.570e-01, 1.200e-02, area=1.680e-02, Iyy=4.160e-04, Izz=4.160e-04, J=8.310e-04)
        S['457.0']['12.5'] = CHS(4.570e-01, 1.250e-02, area=1.750e-02, Iyy=4.310e-04, Izz=4.310e-04, J=8.630e-04)
        S['457.0']['16.0'] = CHS(4.570e-01, 1.600e-02, area=2.220e-02, Iyy=5.400e-04, Izz=5.400e-04, J=1.080e-03)
        S['508.0'] = {}
        S['508.0']['6.0'] = CHS(5.080e-01, 6.000e-03, area=9.460e-03, Iyy=2.980e-04, Izz=2.980e-04, J=5.960e-04)
        S['508.0']['6.3'] = CHS(5.080e-01, 6.300e-03, area=9.930e-03, Iyy=3.120e-04, Izz=3.120e-04, J=6.250e-04)
        S['508.0']['8.0'] = CHS(5.080e-01, 8.000e-03, area=1.260e-02, Iyy=3.930e-04, Izz=3.930e-04, J=7.860e-04)
        S['508.0']['10.0'] = CHS(5.080e-01, 1.000e-02, area=1.560e-02, Iyy=4.850e-04, Izz=4.850e-04, J=9.700e-04)
        S['508.0']['12.0'] = CHS(5.080e-01, 1.200e-02, area=1.870e-02, Iyy=5.750e-04, Izz=5.750e-04, J=1.150e-03)
        S['508.0']['12.5'] = CHS(5.080e-01, 1.250e-02, area=1.950e-02, Iyy=5.980e-04, Izz=5.980e-04, J=1.200e-03)
        S['508.0']['16.0'] = CHS(5.080e-01, 1.600e-02, area=2.470e-02, Iyy=7.490e-04, Izz=7.490e-04, J=1.500e-03)

        return S # No. data rows = 106

    @staticmethod
    def CHS_ColdFormed():
        if BlueBook.__CHS_ColdFormed is None:
            BlueBook.__CHS_ColdFormed = BlueBook.__create_CHS_ColdFormed()
        return BlueBook.__CHS_ColdFormed

    __CHS_HotFinished = None

    @staticmethod
    def __create_CHS_HotFinished():
        S = {}
        S['21.3'] = {}
        S['21.3']['2.6'] = CHS(2.130e-02, 2.600e-03, area=1.530e-04, Iyy=6.810e-09, Izz=6.810e-09, J=1.360e-08)
        S['21.3']['2.9'] = CHS(2.130e-02, 2.900e-03, area=1.680e-04, Iyy=7.270e-09, Izz=7.270e-09, J=1.450e-08)
        S['21.3']['3.2'] = CHS(2.130e-02, 3.200e-03, area=1.820e-04, Iyy=7.680e-09, Izz=7.680e-09, J=1.540e-08)
        S['26.9'] = {}
        S['26.9']['2.6'] = CHS(2.690e-02, 2.600e-03, area=1.980e-04, Iyy=1.480e-08, Izz=1.480e-08, J=2.960e-08)
        S['26.9']['2.9'] = CHS(2.690e-02, 2.900e-03, area=2.190e-04, Iyy=1.600e-08, Izz=1.600e-08, J=3.190e-08)
        S['26.9']['3.2'] = CHS(2.690e-02, 3.200e-03, area=2.380e-04, Iyy=1.700e-08, Izz=1.700e-08, J=3.410e-08)
        S['26.9']['3.6'] = CHS(2.690e-02, 3.600e-03, area=2.640e-04, Iyy=1.830e-08, Izz=1.830e-08, J=3.660e-08)
        S['33.7'] = {}
        S['33.7']['2.6'] = CHS(3.370e-02, 2.600e-03, area=2.540e-04, Iyy=3.090e-08, Izz=3.090e-08, J=6.190e-08)
        S['33.7']['2.9'] = CHS(3.370e-02, 2.900e-03, area=2.810e-04, Iyy=3.360e-08, Izz=3.360e-08, J=6.710e-08)
        S['33.7']['3.2'] = CHS(3.370e-02, 3.200e-03, area=3.070e-04, Iyy=3.600e-08, Izz=3.600e-08, J=7.210e-08)
        S['33.7']['3.6'] = CHS(3.370e-02, 3.600e-03, area=3.400e-04, Iyy=3.910e-08, Izz=3.910e-08, J=7.820e-08)
        S['33.7']['4.0'] = CHS(3.370e-02, 4.000e-03, area=3.730e-04, Iyy=4.190e-08, Izz=4.190e-08, J=8.380e-08)
        S['33.7']['4.5'] = CHS(3.370e-02, 4.500e-03, area=4.130e-04, Iyy=4.500e-08, Izz=4.500e-08, J=9.010e-08)
        S['42.4'] = {}
        S['42.4']['2.6'] = CHS(4.240e-02, 2.600e-03, area=3.250e-04, Iyy=6.460e-08, Izz=6.460e-08, J=1.290e-07)
        S['42.4']['2.9'] = CHS(4.240e-02, 2.900e-03, area=3.600e-04, Iyy=7.060e-08, Izz=7.060e-08, J=1.410e-07)
        S['42.4']['3.2'] = CHS(4.240e-02, 3.200e-03, area=3.940e-04, Iyy=7.620e-08, Izz=7.620e-08, J=1.520e-07)
        S['42.4']['3.6'] = CHS(4.240e-02, 3.600e-03, area=4.390e-04, Iyy=8.330e-08, Izz=8.330e-08, J=1.670e-07)
        S['42.4']['4.0'] = CHS(4.240e-02, 4.000e-03, area=4.830e-04, Iyy=8.990e-08, Izz=8.990e-08, J=1.800e-07)
        S['42.4']['4.5'] = CHS(4.240e-02, 4.500e-03, area=5.360e-04, Iyy=9.760e-08, Izz=9.760e-08, J=1.950e-07)
        S['48.3'] = {}
        S['48.3']['2.6'] = CHS(4.830e-02, 2.600e-03, area=3.730e-04, Iyy=9.780e-08, Izz=9.780e-08, J=1.960e-07)
        S['48.3']['2.9'] = CHS(4.830e-02, 2.900e-03, area=4.140e-04, Iyy=1.070e-07, Izz=1.070e-07, J=2.140e-07)
        S['48.3']['3.2'] = CHS(4.830e-02, 3.200e-03, area=4.530e-04, Iyy=1.160e-07, Izz=1.160e-07, J=2.320e-07)
        S['48.3']['3.6'] = CHS(4.830e-02, 3.600e-03, area=5.060e-04, Iyy=1.270e-07, Izz=1.270e-07, J=2.540e-07)
        S['48.3']['4.0'] = CHS(4.830e-02, 4.000e-03, area=5.570e-04, Iyy=1.380e-07, Izz=1.380e-07, J=2.750e-07)
        S['48.3']['4.5'] = CHS(4.830e-02, 4.500e-03, area=6.190e-04, Iyy=1.500e-07, Izz=1.500e-07, J=3.000e-07)
        S['48.3']['5.0'] = CHS(4.830e-02, 5.000e-03, area=6.800e-04, Iyy=1.620e-07, Izz=1.620e-07, J=3.230e-07)
        S['48.3']['5.6'] = CHS(4.830e-02, 5.600e-03, area=7.510e-04, Iyy=1.740e-07, Izz=1.740e-07, J=3.480e-07)
        S['48.3']['6.3'] = CHS(4.830e-02, 6.300e-03, area=8.310e-04, Iyy=1.870e-07, Izz=1.870e-07, J=3.750e-07)
        S['60.3'] = {}
        S['60.3']['2.6'] = CHS(6.030e-02, 2.600e-03, area=4.710e-04, Iyy=1.970e-07, Izz=1.970e-07, J=3.930e-07)
        S['60.3']['2.9'] = CHS(6.030e-02, 2.900e-03, area=5.230e-04, Iyy=2.160e-07, Izz=2.160e-07, J=4.320e-07)
        S['60.3']['3.2'] = CHS(6.030e-02, 3.200e-03, area=5.740e-04, Iyy=2.350e-07, Izz=2.350e-07, J=4.690e-07)
        S['60.3']['3.6'] = CHS(6.030e-02, 3.600e-03, area=6.410e-04, Iyy=2.590e-07, Izz=2.590e-07, J=5.170e-07)
        S['60.3']['4.0'] = CHS(6.030e-02, 4.000e-03, area=7.070e-04, Iyy=2.820e-07, Izz=2.820e-07, J=5.630e-07)
        S['60.3']['4.5'] = CHS(6.030e-02, 4.500e-03, area=7.890e-04, Iyy=3.090e-07, Izz=3.090e-07, J=6.180e-07)
        S['60.3']['5.0'] = CHS(6.030e-02, 5.000e-03, area=8.690e-04, Iyy=3.350e-07, Izz=3.350e-07, J=6.700e-07)
        S['60.3']['5.6'] = CHS(6.030e-02, 5.600e-03, area=9.620e-04, Iyy=3.640e-07, Izz=3.640e-07, J=7.270e-07)
        S['60.3']['6.3'] = CHS(6.030e-02, 6.300e-03, area=1.070e-03, Iyy=3.950e-07, Izz=3.950e-07, J=7.900e-07)
        S['60.3']['7.1'] = CHS(6.030e-02, 7.100e-03, area=1.190e-03, Iyy=4.270e-07, Izz=4.270e-07, J=8.550e-07)
        S['60.3']['8.0'] = CHS(6.030e-02, 8.000e-03, area=1.310e-03, Iyy=4.600e-07, Izz=4.600e-07, J=9.200e-07)
        S['76.1'] = {}
        S['76.1']['2.9'] = CHS(7.610e-02, 2.900e-03, area=6.670e-04, Iyy=4.470e-07, Izz=4.470e-07, J=8.950e-07)
        S['76.1']['3.2'] = CHS(7.610e-02, 3.200e-03, area=7.330e-04, Iyy=4.880e-07, Izz=4.880e-07, J=9.760e-07)
        S['76.1']['3.6'] = CHS(7.610e-02, 3.600e-03, area=8.200e-04, Iyy=5.400e-07, Izz=5.400e-07, J=1.080e-06)
        S['76.1']['4.0'] = CHS(7.610e-02, 4.000e-03, area=9.060e-04, Iyy=5.910e-07, Izz=5.910e-07, J=1.180e-06)
        S['76.1']['4.5'] = CHS(7.610e-02, 4.500e-03, area=1.010e-03, Iyy=6.510e-07, Izz=6.510e-07, J=1.300e-06)
        S['76.1']['5.0'] = CHS(7.610e-02, 5.000e-03, area=1.120e-03, Iyy=7.090e-07, Izz=7.090e-07, J=1.420e-06)
        S['76.1']['5.6'] = CHS(7.610e-02, 5.600e-03, area=1.240e-03, Iyy=7.750e-07, Izz=7.750e-07, J=1.550e-06)
        S['76.1']['6.3'] = CHS(7.610e-02, 6.300e-03, area=1.380e-03, Iyy=8.480e-07, Izz=8.480e-07, J=1.700e-06)
        S['76.1']['7.1'] = CHS(7.610e-02, 7.100e-03, area=1.540e-03, Iyy=9.260e-07, Izz=9.260e-07, J=1.850e-06)
        S['76.1']['8.0'] = CHS(7.610e-02, 8.000e-03, area=1.710e-03, Iyy=1.010e-06, Izz=1.010e-06, J=2.010e-06)
        S['88.9'] = {}
        S['88.9']['3.2'] = CHS(8.890e-02, 3.200e-03, area=8.620e-04, Iyy=7.920e-07, Izz=7.920e-07, J=1.580e-06)
        S['88.9']['3.6'] = CHS(8.890e-02, 3.600e-03, area=9.650e-04, Iyy=8.790e-07, Izz=8.790e-07, J=1.760e-06)
        S['88.9']['4.0'] = CHS(8.890e-02, 4.000e-03, area=1.070e-03, Iyy=9.630e-07, Izz=9.630e-07, J=1.930e-06)
        S['88.9']['4.5'] = CHS(8.890e-02, 4.500e-03, area=1.190e-03, Iyy=1.070e-06, Izz=1.070e-06, J=2.130e-06)
        S['88.9']['5.0'] = CHS(8.890e-02, 5.000e-03, area=1.320e-03, Iyy=1.160e-06, Izz=1.160e-06, J=2.330e-06)
        S['88.9']['5.6'] = CHS(8.890e-02, 5.600e-03, area=1.470e-03, Iyy=1.280e-06, Izz=1.280e-06, J=2.550e-06)
        S['88.9']['6.3'] = CHS(8.890e-02, 6.300e-03, area=1.630e-03, Iyy=1.400e-06, Izz=1.400e-06, J=2.800e-06)
        S['88.9']['7.1'] = CHS(8.890e-02, 7.100e-03, area=1.820e-03, Iyy=1.540e-06, Izz=1.540e-06, J=3.080e-06)
        S['88.9']['8.0'] = CHS(8.890e-02, 8.000e-03, area=2.030e-03, Iyy=1.680e-06, Izz=1.680e-06, J=3.360e-06)
        S['88.9']['10.0'] = CHS(8.890e-02, 1.000e-02, area=2.480e-03, Iyy=1.960e-06, Izz=1.960e-06, J=3.920e-06)
        S['101.6'] = {}
        S['101.6']['3.6'] = CHS(1.016e-01, 3.600e-03, area=1.110e-03, Iyy=1.330e-06, Izz=1.330e-06, J=2.660e-06)
        S['101.6']['4.0'] = CHS(1.016e-01, 4.000e-03, area=1.230e-03, Iyy=1.460e-06, Izz=1.460e-06, J=2.930e-06)
        S['101.6']['4.5'] = CHS(1.016e-01, 4.500e-03, area=1.370e-03, Iyy=1.620e-06, Izz=1.620e-06, J=3.240e-06)
        S['101.6']['5.0'] = CHS(1.016e-01, 5.000e-03, area=1.520e-03, Iyy=1.770e-06, Izz=1.770e-06, J=3.550e-06)
        S['101.6']['5.6'] = CHS(1.016e-01, 5.600e-03, area=1.690e-03, Iyy=1.950e-06, Izz=1.950e-06, J=3.900e-06)
        S['101.6']['6.3'] = CHS(1.016e-01, 6.300e-03, area=1.890e-03, Iyy=2.150e-06, Izz=2.150e-06, J=4.300e-06)
        S['101.6']['7.1'] = CHS(1.016e-01, 7.100e-03, area=2.110e-03, Iyy=2.370e-06, Izz=2.370e-06, J=4.730e-06)
        S['101.6']['8.0'] = CHS(1.016e-01, 8.000e-03, area=2.350e-03, Iyy=2.600e-06, Izz=2.600e-06, J=5.190e-06)
        S['101.6']['10.0'] = CHS(1.016e-01, 1.000e-02, area=2.880e-03, Iyy=3.050e-06, Izz=3.050e-06, J=6.110e-06)
        S['114.3'] = {}
        S['114.3']['3.6'] = CHS(1.143e-01, 3.600e-03, area=1.250e-03, Iyy=1.920e-06, Izz=1.920e-06, J=3.840e-06)
        S['114.3']['4.0'] = CHS(1.143e-01, 4.000e-03, area=1.390e-03, Iyy=2.110e-06, Izz=2.110e-06, J=4.220e-06)
        S['114.3']['4.5'] = CHS(1.143e-01, 4.500e-03, area=1.550e-03, Iyy=2.340e-06, Izz=2.340e-06, J=4.690e-06)
        S['114.3']['5.0'] = CHS(1.143e-01, 5.000e-03, area=1.720e-03, Iyy=2.570e-06, Izz=2.570e-06, J=5.140e-06)
        S['114.3']['5.6'] = CHS(1.143e-01, 5.600e-03, area=1.910e-03, Iyy=2.830e-06, Izz=2.830e-06, J=5.660e-06)
        S['114.3']['6.3'] = CHS(1.143e-01, 6.300e-03, area=2.140e-03, Iyy=3.130e-06, Izz=3.130e-06, J=6.250e-06)
        S['114.3']['7.1'] = CHS(1.143e-01, 7.100e-03, area=2.390e-03, Iyy=3.450e-06, Izz=3.450e-06, J=6.900e-06)
        S['114.3']['8.0'] = CHS(1.143e-01, 8.000e-03, area=2.670e-03, Iyy=3.790e-06, Izz=3.790e-06, J=7.590e-06)
        S['114.3']['10.0'] = CHS(1.143e-01, 1.000e-02, area=3.280e-03, Iyy=4.500e-06, Izz=4.500e-06, J=8.990e-06)
        S['139.7'] = {}
        S['139.7']['3.6'] = CHS(1.397e-01, 3.600e-03, area=1.540e-03, Iyy=3.570e-06, Izz=3.570e-06, J=7.130e-06)
        S['139.7']['4.0'] = CHS(1.397e-01, 4.000e-03, area=1.710e-03, Iyy=3.930e-06, Izz=3.930e-06, J=7.860e-06)
        S['139.7']['4.5'] = CHS(1.397e-01, 4.500e-03, area=1.910e-03, Iyy=4.370e-06, Izz=4.370e-06, J=8.740e-06)
        S['139.7']['5.0'] = CHS(1.397e-01, 5.000e-03, area=2.120e-03, Iyy=4.810e-06, Izz=4.810e-06, J=9.610e-06)
        S['139.7']['5.6'] = CHS(1.397e-01, 5.600e-03, area=2.360e-03, Iyy=5.310e-06, Izz=5.310e-06, J=1.060e-05)
        S['139.7']['6.3'] = CHS(1.397e-01, 6.300e-03, area=2.640e-03, Iyy=5.890e-06, Izz=5.890e-06, J=1.180e-05)
        S['139.7']['7.1'] = CHS(1.397e-01, 7.100e-03, area=2.960e-03, Iyy=6.520e-06, Izz=6.520e-06, J=1.300e-05)
        S['139.7']['8.0'] = CHS(1.397e-01, 8.000e-03, area=3.310e-03, Iyy=7.200e-06, Izz=7.200e-06, J=1.440e-05)
        S['139.7']['10.0'] = CHS(1.397e-01, 1.000e-02, area=4.070e-03, Iyy=8.620e-06, Izz=8.620e-06, J=1.720e-05)
        S['168.3'] = {}
        S['168.3']['5.0'] = CHS(1.683e-01, 5.000e-03, area=2.570e-03, Iyy=8.560e-06, Izz=8.560e-06, J=1.710e-05)
        S['168.3']['5.6'] = CHS(1.683e-01, 5.600e-03, area=2.860e-03, Iyy=9.480e-06, Izz=9.480e-06, J=1.900e-05)
        S['168.3']['6.3'] = CHS(1.683e-01, 6.300e-03, area=3.210e-03, Iyy=1.050e-05, Izz=1.050e-05, J=2.110e-05)
        S['168.3']['7.1'] = CHS(1.683e-01, 7.100e-03, area=3.600e-03, Iyy=1.170e-05, Izz=1.170e-05, J=2.340e-05)
        S['168.3']['8.0'] = CHS(1.683e-01, 8.000e-03, area=4.030e-03, Iyy=1.300e-05, Izz=1.300e-05, J=2.600e-05)
        S['168.3']['10.0'] = CHS(1.683e-01, 1.000e-02, area=4.970e-03, Iyy=1.560e-05, Izz=1.560e-05, J=3.130e-05)
        S['168.3']['11.0'] = CHS(1.683e-01, 1.100e-02, area=5.440e-03, Iyy=1.690e-05, Izz=1.690e-05, J=3.380e-05)
        S['168.3']['12.5'] = CHS(1.683e-01, 1.250e-02, area=6.120e-03, Iyy=1.870e-05, Izz=1.870e-05, J=3.740e-05)
        S['193.7'] = {}
        S['193.7']['5.0'] = CHS(1.937e-01, 5.000e-03, area=2.960e-03, Iyy=1.320e-05, Izz=1.320e-05, J=2.640e-05)
        S['193.7']['5.6'] = CHS(1.937e-01, 5.600e-03, area=3.310e-03, Iyy=1.460e-05, Izz=1.460e-05, J=2.930e-05)
        S['193.7']['6.3'] = CHS(1.937e-01, 6.300e-03, area=3.710e-03, Iyy=1.630e-05, Izz=1.630e-05, J=3.260e-05)
        S['193.7']['7.1'] = CHS(1.937e-01, 7.100e-03, area=4.160e-03, Iyy=1.810e-05, Izz=1.810e-05, J=3.630e-05)
        S['193.7']['8.0'] = CHS(1.937e-01, 8.000e-03, area=4.670e-03, Iyy=2.020e-05, Izz=2.020e-05, J=4.030e-05)
        S['193.7']['10.0'] = CHS(1.937e-01, 1.000e-02, area=5.770e-03, Iyy=2.440e-05, Izz=2.440e-05, J=4.880e-05)
        S['193.7']['11.0'] = CHS(1.937e-01, 1.100e-02, area=6.310e-03, Iyy=2.640e-05, Izz=2.640e-05, J=5.290e-05)
        S['193.7']['12.5'] = CHS(1.937e-01, 1.250e-02, area=7.120e-03, Iyy=2.930e-05, Izz=2.930e-05, J=5.870e-05)
        S['193.7']['14.2'] = CHS(1.937e-01, 1.420e-02, area=8.010e-03, Iyy=3.240e-05, Izz=3.240e-05, J=6.490e-05)
        S['193.7']['16.0'] = CHS(1.937e-01, 1.600e-02, area=8.930e-03, Iyy=3.550e-05, Izz=3.550e-05, J=7.110e-05)
        S['219.1'] = {}
        S['219.1']['4.5'] = CHS(2.191e-01, 4.500e-03, area=3.030e-03, Iyy=1.750e-05, Izz=1.750e-05, J=3.490e-05)
        S['219.1']['5.0'] = CHS(2.191e-01, 5.000e-03, area=3.360e-03, Iyy=1.930e-05, Izz=1.930e-05, J=3.860e-05)
        S['219.1']['5.6'] = CHS(2.191e-01, 5.600e-03, area=3.760e-03, Iyy=2.140e-05, Izz=2.140e-05, J=4.280e-05)
        S['219.1']['6.3'] = CHS(2.191e-01, 6.300e-03, area=4.210e-03, Iyy=2.390e-05, Izz=2.390e-05, J=4.770e-05)
        S['219.1']['7.1'] = CHS(2.191e-01, 7.100e-03, area=4.730e-03, Iyy=2.660e-05, Izz=2.660e-05, J=5.320e-05)
        S['219.1']['8.0'] = CHS(2.191e-01, 8.000e-03, area=5.310e-03, Iyy=2.960e-05, Izz=2.960e-05, J=5.920e-05)
        S['219.1']['10.0'] = CHS(2.191e-01, 1.000e-02, area=6.570e-03, Iyy=3.600e-05, Izz=3.600e-05, J=7.200e-05)
        S['219.1']['11.0'] = CHS(2.191e-01, 1.100e-02, area=7.190e-03, Iyy=3.900e-05, Izz=3.900e-05, J=7.810e-05)
        S['219.1']['12.5'] = CHS(2.191e-01, 1.250e-02, area=8.110e-03, Iyy=4.340e-05, Izz=4.340e-05, J=8.690e-05)
        S['219.1']['14.2'] = CHS(2.191e-01, 1.420e-02, area=9.140e-03, Iyy=4.820e-05, Izz=4.820e-05, J=9.640e-05)
        S['219.1']['16.0'] = CHS(2.191e-01, 1.600e-02, area=1.020e-02, Iyy=5.300e-05, Izz=5.300e-05, J=1.060e-04)
        S['244.5'] = {}
        S['244.5']['5.0'] = CHS(2.445e-01, 5.000e-03, area=3.760e-03, Iyy=2.700e-05, Izz=2.700e-05, J=5.400e-05)
        S['244.5']['5.6'] = CHS(2.445e-01, 5.600e-03, area=4.200e-03, Iyy=3.000e-05, Izz=3.000e-05, J=6.000e-05)
        S['244.5']['6.3'] = CHS(2.445e-01, 6.300e-03, area=4.710e-03, Iyy=3.350e-05, Izz=3.350e-05, J=6.690e-05)
        S['244.5']['7.1'] = CHS(2.445e-01, 7.100e-03, area=5.300e-03, Iyy=3.730e-05, Izz=3.730e-05, J=7.470e-05)
        S['244.5']['8.0'] = CHS(2.445e-01, 8.000e-03, area=5.940e-03, Iyy=4.160e-05, Izz=4.160e-05, J=8.320e-05)
        S['244.5']['10.0'] = CHS(2.445e-01, 1.000e-02, area=7.370e-03, Iyy=5.070e-05, Izz=5.070e-05, J=1.010e-04)
        S['244.5']['11.0'] = CHS(2.445e-01, 1.100e-02, area=8.070e-03, Iyy=5.510e-05, Izz=5.510e-05, J=1.100e-04)
        S['244.5']['12.5'] = CHS(2.445e-01, 1.250e-02, area=9.110e-03, Iyy=6.150e-05, Izz=6.150e-05, J=1.230e-04)
        S['244.5']['14.2'] = CHS(2.445e-01, 1.420e-02, area=1.030e-02, Iyy=6.840e-05, Izz=6.840e-05, J=1.370e-04)
        S['244.5']['16.0'] = CHS(2.445e-01, 1.600e-02, area=1.150e-02, Iyy=7.530e-05, Izz=7.530e-05, J=1.510e-04)
        S['273.0'] = {}
        S['273.0']['5.0'] = CHS(2.730e-01, 5.000e-03, area=4.210e-03, Iyy=3.780e-05, Izz=3.780e-05, J=7.560e-05)
        S['273.0']['5.6'] = CHS(2.730e-01, 5.600e-03, area=4.700e-03, Iyy=4.210e-05, Izz=4.210e-05, J=8.410e-05)
        S['273.0']['6.3'] = CHS(2.730e-01, 6.300e-03, area=5.280e-03, Iyy=4.700e-05, Izz=4.700e-05, J=9.390e-05)
        S['273.0']['7.1'] = CHS(2.730e-01, 7.100e-03, area=5.930e-03, Iyy=5.240e-05, Izz=5.240e-05, J=1.050e-04)
        S['273.0']['8.0'] = CHS(2.730e-01, 8.000e-03, area=6.660e-03, Iyy=5.850e-05, Izz=5.850e-05, J=1.170e-04)
        S['273.0']['10.0'] = CHS(2.730e-01, 1.000e-02, area=8.260e-03, Iyy=7.150e-05, Izz=7.150e-05, J=1.430e-04)
        S['273.0']['11.0'] = CHS(2.730e-01, 1.100e-02, area=9.050e-03, Iyy=7.780e-05, Izz=7.780e-05, J=1.560e-04)
        S['273.0']['12.5'] = CHS(2.730e-01, 1.250e-02, area=1.020e-02, Iyy=8.700e-05, Izz=8.700e-05, J=1.740e-04)
        S['273.0']['14.2'] = CHS(2.730e-01, 1.420e-02, area=1.150e-02, Iyy=9.700e-05, Izz=9.700e-05, J=1.940e-04)
        S['273.0']['16.0'] = CHS(2.730e-01, 1.600e-02, area=1.290e-02, Iyy=1.070e-04, Izz=1.070e-04, J=2.140e-04)
        S['323.9'] = {}
        S['323.9']['5.0'] = CHS(3.239e-01, 5.000e-03, area=5.010e-03, Iyy=6.370e-05, Izz=6.370e-05, J=1.270e-04)
        S['323.9']['5.6'] = CHS(3.239e-01, 5.600e-03, area=5.600e-03, Iyy=7.090e-05, Izz=7.090e-05, J=1.420e-04)
        S['323.9']['6.3'] = CHS(3.239e-01, 6.300e-03, area=6.290e-03, Iyy=7.930e-05, Izz=7.930e-05, J=1.590e-04)
        S['323.9']['7.1'] = CHS(3.239e-01, 7.100e-03, area=7.070e-03, Iyy=8.870e-05, Izz=8.870e-05, J=1.770e-04)
        S['323.9']['8.0'] = CHS(3.239e-01, 8.000e-03, area=7.940e-03, Iyy=9.910e-05, Izz=9.910e-05, J=1.980e-04)
        S['323.9']['10.0'] = CHS(3.239e-01, 1.000e-02, area=9.860e-03, Iyy=1.220e-04, Izz=1.220e-04, J=2.430e-04)
        S['323.9']['11.0'] = CHS(3.239e-01, 1.100e-02, area=1.080e-02, Iyy=1.320e-04, Izz=1.320e-04, J=2.650e-04)
        S['323.9']['12.5'] = CHS(3.239e-01, 1.250e-02, area=1.220e-02, Iyy=1.480e-04, Izz=1.480e-04, J=2.970e-04)
        S['323.9']['14.2'] = CHS(3.239e-01, 1.420e-02, area=1.380e-02, Iyy=1.660e-04, Izz=1.660e-04, J=3.320e-04)
        S['323.9']['16.0'] = CHS(3.239e-01, 1.600e-02, area=1.550e-02, Iyy=1.840e-04, Izz=1.840e-04, J=3.680e-04)
        S['355.6'] = {}
        S['355.6']['6.3'] = CHS(3.556e-01, 6.300e-03, area=6.910e-03, Iyy=1.050e-04, Izz=1.050e-04, J=2.110e-04)
        S['355.6']['7.1'] = CHS(3.556e-01, 7.100e-03, area=7.770e-03, Iyy=1.180e-04, Izz=1.180e-04, J=2.360e-04)
        S['355.6']['8.0'] = CHS(3.556e-01, 8.000e-03, area=8.740e-03, Iyy=1.320e-04, Izz=1.320e-04, J=2.640e-04)
        S['355.6']['10.0'] = CHS(3.556e-01, 1.000e-02, area=1.090e-02, Iyy=1.620e-04, Izz=1.620e-04, J=3.240e-04)
        S['355.6']['11.0'] = CHS(3.556e-01, 1.100e-02, area=1.190e-02, Iyy=1.770e-04, Izz=1.770e-04, J=3.540e-04)
        S['355.6']['12.5'] = CHS(3.556e-01, 1.250e-02, area=1.350e-02, Iyy=1.990e-04, Izz=1.990e-04, J=3.970e-04)
        S['355.6']['14.2'] = CHS(3.556e-01, 1.420e-02, area=1.520e-02, Iyy=2.220e-04, Izz=2.220e-04, J=4.450e-04)
        S['355.6']['16.0'] = CHS(3.556e-01, 1.600e-02, area=1.710e-02, Iyy=2.470e-04, Izz=2.470e-04, J=4.930e-04)
        S['406.4'] = {}
        S['406.4']['10.0'] = CHS(4.064e-01, 1.000e-02, area=1.250e-02, Iyy=2.450e-04, Izz=2.450e-04, J=4.900e-04)
        S['406.4']['11.0'] = CHS(4.064e-01, 1.100e-02, area=1.370e-02, Iyy=2.670e-04, Izz=2.670e-04, J=5.340e-04)
        S['406.4']['12.5'] = CHS(4.064e-01, 1.250e-02, area=1.550e-02, Iyy=3.000e-04, Izz=3.000e-04, J=6.010e-04)
        S['406.4']['14.2'] = CHS(4.064e-01, 1.420e-02, area=1.750e-02, Iyy=3.370e-04, Izz=3.370e-04, J=6.740e-04)
        S['406.4']['16.0'] = CHS(4.064e-01, 1.600e-02, area=1.960e-02, Iyy=3.740e-04, Izz=3.740e-04, J=7.490e-04)
        S['457.0'] = {}
        S['457.0']['10.0'] = CHS(4.570e-01, 1.000e-02, area=1.400e-02, Iyy=3.510e-04, Izz=3.510e-04, J=7.020e-04)
        S['457.0']['11.0'] = CHS(4.570e-01, 1.100e-02, area=1.540e-02, Iyy=3.830e-04, Izz=3.830e-04, J=7.670e-04)
        S['457.0']['12.5'] = CHS(4.570e-01, 1.250e-02, area=1.750e-02, Iyy=4.310e-04, Izz=4.310e-04, J=8.630e-04)
        S['457.0']['14.2'] = CHS(4.570e-01, 1.420e-02, area=1.980e-02, Iyy=4.850e-04, Izz=4.850e-04, J=9.690e-04)
        S['457.0']['16.0'] = CHS(4.570e-01, 1.600e-02, area=2.220e-02, Iyy=5.400e-04, Izz=5.400e-04, J=1.080e-03)
        S['508.0'] = {}
        S['508.0']['10.0'] = CHS(5.080e-01, 1.000e-02, area=1.560e-02, Iyy=4.850e-04, Izz=4.850e-04, J=9.700e-04)
        S['508.0']['11.0'] = CHS(5.080e-01, 1.100e-02, area=1.720e-02, Iyy=5.310e-04, Izz=5.310e-04, J=1.060e-03)
        S['508.0']['12.5'] = CHS(5.080e-01, 1.250e-02, area=1.950e-02, Iyy=5.980e-04, Izz=5.980e-04, J=1.200e-03)
        S['508.0']['14.2'] = CHS(5.080e-01, 1.420e-02, area=2.200e-02, Iyy=6.720e-04, Izz=6.720e-04, J=1.340e-03)
        S['508.0']['16.0'] = CHS(5.080e-01, 1.600e-02, area=2.470e-02, Iyy=7.490e-04, Izz=7.490e-04, J=1.500e-03)

        return S # No. data rows = 168

    @staticmethod
    def CHS_HotFinished():
        if BlueBook.__CHS_HotFinished is None:
            BlueBook.__CHS_HotFinished = BlueBook.__create_CHS_HotFinished()
        return BlueBook.__CHS_HotFinished

    __RHS_ColdFormed = None

    @staticmethod
    def __create_RHS_ColdFormed():
        S = {}
        S['50  x  25'] = {}
        S['50  x  25']['2.0'] = RHS(2.500e-02, 5.000e-02, 2.000e-03, 5.000e-03, area=2.740e-04, Iyy=8.380e-08, Izz=2.810e-08, J=7.060e-08)
        S['50  x  25']['2.5'] = RHS(2.500e-02, 5.000e-02, 2.500e-03, 6.250e-03, area=3.340e-04, Iyy=9.890e-08, Izz=3.280e-08, J=8.430e-08)
        S['50  x  25']['3.0'] = RHS(2.500e-02, 5.000e-02, 3.000e-03, 7.500e-03, area=3.910e-04, Iyy=1.120e-07, Izz=3.670e-08, J=9.640e-08)
        S['50  x  30'] = {}
        S['50  x  30']['2.0'] = RHS(3.000e-02, 5.000e-02, 2.000e-03, 5.000e-03, area=2.940e-04, Iyy=9.540e-08, Izz=4.290e-08, J=9.770e-08)
        S['50  x  30']['2.5'] = RHS(3.000e-02, 5.000e-02, 2.500e-03, 6.250e-03, area=3.590e-04, Iyy=1.130e-07, Izz=5.050e-08, J=1.170e-07)
        S['50  x  30']['3.0'] = RHS(3.000e-02, 5.000e-02, 3.000e-03, 7.500e-03, area=4.210e-04, Iyy=1.280e-07, Izz=5.700e-08, J=1.350e-07)
        S['50  x  30']['4.0'] = RHS(3.000e-02, 5.000e-02, 4.000e-03, 1.000e-02, area=5.350e-04, Iyy=1.530e-07, Izz=6.690e-08, J=1.650e-07)
        S['60  x  40'] = {}
        S['60  x  40']['2.5'] = RHS(4.000e-02, 6.000e-02, 2.500e-03, 6.250e-03, area=4.590e-04, Iyy=2.210e-07, Izz=1.170e-07, J=2.510e-07)
        S['60  x  40']['3.0'] = RHS(4.000e-02, 6.000e-02, 3.000e-03, 7.500e-03, area=5.410e-04, Iyy=2.540e-07, Izz=1.340e-07, J=2.930e-07)
        S['60  x  40']['4.0'] = RHS(4.000e-02, 6.000e-02, 4.000e-03, 1.000e-02, area=6.950e-04, Iyy=3.100e-07, Izz=1.630e-07, J=3.670e-07)
        S['60  x  40']['5.0'] = RHS(4.000e-02, 6.000e-02, 5.000e-03, 1.250e-02, area=8.360e-04, Iyy=3.530e-07, Izz=1.840e-07, J=4.280e-07)
        S['70  x  40'] = {}
        S['70  x  40']['3.0'] = RHS(4.000e-02, 7.000e-02, 3.000e-03, 7.500e-03, area=6.010e-04, Iyy=3.730e-07, Izz=1.550e-07, J=3.650e-07)
        S['70  x  40']['4.0'] = RHS(4.000e-02, 7.000e-02, 4.000e-03, 1.000e-02, area=7.750e-04, Iyy=4.600e-07, Izz=1.890e-07, J=4.580e-07)
        S['70  x  40']['5.0'] = RHS(4.000e-02, 7.000e-02, 5.000e-03, 1.250e-02, area=9.360e-04, Iyy=5.290e-07, Izz=2.150e-07, J=5.380e-07)
        S['70  x  50'] = {}
        S['70  x  50']['3.0'] = RHS(5.000e-02, 7.000e-02, 3.000e-03, 7.500e-03, area=6.610e-04, Iyy=4.410e-07, Izz=2.610e-07, J=5.360e-07)
        S['70  x  50']['4.0'] = RHS(5.000e-02, 7.000e-02, 4.000e-03, 1.000e-02, area=8.550e-04, Iyy=5.470e-07, Izz=3.220e-07, J=6.810e-07)
        S['70  x  50']['5.0'] = RHS(5.000e-02, 7.000e-02, 5.000e-03, 1.250e-02, area=1.040e-03, Iyy=6.350e-07, Izz=3.720e-07, J=8.080e-07)
        S['80  x  40'] = {}
        S['80  x  40']['3.0'] = RHS(4.000e-02, 8.000e-02, 3.000e-03, 7.500e-03, area=6.610e-04, Iyy=5.230e-07, Izz=1.760e-07, J=4.390e-07)
        S['80  x  40']['4.0'] = RHS(4.000e-02, 8.000e-02, 4.000e-03, 1.000e-02, area=8.550e-04, Iyy=6.480e-07, Izz=2.150e-07, J=5.520e-07)
        S['80  x  40']['5.0'] = RHS(4.000e-02, 8.000e-02, 5.000e-03, 1.250e-02, area=1.040e-03, Iyy=7.510e-07, Izz=2.460e-07, J=6.500e-07)
        S['80  x  50'] = {}
        S['80  x  50']['3.0'] = RHS(5.000e-02, 8.000e-02, 3.000e-03, 7.500e-03, area=7.210e-04, Iyy=6.110e-07, Izz=2.940e-07, J=6.500e-07)
        S['80  x  50']['4.0'] = RHS(5.000e-02, 8.000e-02, 4.000e-03, 1.000e-02, area=9.350e-04, Iyy=7.640e-07, Izz=3.650e-07, J=8.270e-07)
        S['80  x  50']['5.0'] = RHS(5.000e-02, 8.000e-02, 5.000e-03, 1.250e-02, area=1.140e-03, Iyy=8.920e-07, Izz=4.230e-07, J=9.840e-07)
        S['80  x  60'] = {}
        S['80  x  60']['3.0'] = RHS(6.000e-02, 8.000e-02, 3.000e-03, 7.500e-03, area=7.810e-04, Iyy=7.000e-07, Izz=4.490e-07, J=8.830e-07)
        S['80  x  60']['3.5'] = RHS(6.000e-02, 8.000e-02, 3.500e-03, 8.750e-03, area=8.990e-04, Iyy=7.930e-07, Izz=5.070e-07, J=1.010e-06)
        S['80  x  60']['4.0'] = RHS(6.000e-02, 8.000e-02, 4.000e-03, 1.000e-02, area=1.010e-03, Iyy=8.790e-07, Izz=5.610e-07, J=1.130e-06)
        S['80  x  60']['5.0'] = RHS(6.000e-02, 8.000e-02, 5.000e-03, 1.250e-02, area=1.240e-03, Iyy=1.030e-06, Izz=6.570e-07, J=1.360e-06)
        S['90  x  50'] = {}
        S['90  x  50']['3.0'] = RHS(5.000e-02, 9.000e-02, 3.000e-03, 7.500e-03, area=7.810e-04, Iyy=8.190e-07, Izz=3.270e-07, J=7.670e-07)
        S['90  x  50']['4.0'] = RHS(5.000e-02, 9.000e-02, 4.000e-03, 1.000e-02, area=1.010e-03, Iyy=1.030e-06, Izz=4.070e-07, J=9.770e-07)
        S['90  x  50']['5.0'] = RHS(5.000e-02, 9.000e-02, 5.000e-03, 1.250e-02, area=1.240e-03, Iyy=1.210e-06, Izz=4.740e-07, J=1.160e-06)
        S['100  x  40'] = {}
        S['100  x  40']['3.0'] = RHS(4.000e-02, 1.000e-01, 3.000e-03, 7.500e-03, area=7.810e-04, Iyy=9.230e-07, Izz=2.170e-07, J=5.900e-07)
        S['100  x  40']['4.0'] = RHS(4.000e-02, 1.000e-01, 4.000e-03, 1.000e-02, area=1.010e-03, Iyy=1.160e-06, Izz=2.670e-07, J=7.450e-07)
        S['100  x  40']['5.0'] = RHS(4.000e-02, 1.000e-01, 5.000e-03, 1.250e-02, area=1.240e-03, Iyy=1.360e-06, Izz=3.080e-07, J=8.790e-07)
        S['100  x  50'] = {}
        S['100  x  50']['3.0'] = RHS(5.000e-02, 1.000e-01, 3.000e-03, 7.500e-03, area=8.410e-04, Iyy=1.060e-06, Izz=3.610e-07, J=8.860e-07)
        S['100  x  50']['4.0'] = RHS(5.000e-02, 1.000e-01, 4.000e-03, 1.000e-02, area=1.090e-03, Iyy=1.340e-06, Izz=4.490e-07, J=1.130e-06)
        S['100  x  50']['5.0'] = RHS(5.000e-02, 1.000e-01, 5.000e-03, 1.250e-02, area=1.340e-03, Iyy=1.580e-06, Izz=5.250e-07, J=1.350e-06)
        S['100  x  50']['6.0'] = RHS(5.000e-02, 1.000e-01, 6.000e-03, 1.500e-02, area=1.560e-03, Iyy=1.790e-06, Izz=5.870e-07, J=1.540e-06)
        S['100  x  60'] = {}
        S['100  x  60']['3.0'] = RHS(6.000e-02, 1.000e-01, 3.000e-03, 7.500e-03, area=9.010e-04, Iyy=1.210e-06, Izz=5.460e-07, J=1.220e-06)
        S['100  x  60']['3.5'] = RHS(6.000e-02, 1.000e-01, 3.500e-03, 8.750e-03, area=1.040e-03, Iyy=1.370e-06, Izz=6.190e-07, J=1.390e-06)
        S['100  x  60']['4.0'] = RHS(6.000e-02, 1.000e-01, 4.000e-03, 1.000e-02, area=1.170e-03, Iyy=1.530e-06, Izz=6.870e-07, J=1.560e-06)
        S['100  x  60']['5.0'] = RHS(6.000e-02, 1.000e-01, 5.000e-03, 1.250e-02, area=1.440e-03, Iyy=1.810e-06, Izz=8.080e-07, J=1.880e-06)
        S['100  x  60']['6.0'] = RHS(6.000e-02, 1.000e-01, 6.000e-03, 1.500e-02, area=1.680e-03, Iyy=2.050e-06, Izz=9.120e-07, J=2.160e-06)
        S['100  x  80'] = {}
        S['100  x  80']['3.0'] = RHS(8.000e-02, 1.000e-01, 3.000e-03, 7.500e-03, area=1.020e-03, Iyy=1.490e-06, Izz=1.060e-06, J=1.960e-06)
        S['100  x  80']['4.0'] = RHS(8.000e-02, 1.000e-01, 4.000e-03, 1.000e-02, area=1.330e-03, Iyy=1.890e-06, Izz=1.340e-06, J=2.540e-06)
        S['100  x  80']['5.0'] = RHS(8.000e-02, 1.000e-01, 5.000e-03, 1.250e-02, area=1.640e-03, Iyy=2.260e-06, Izz=1.600e-06, J=3.080e-06)
        S['100  x  80']['6.0'] = RHS(8.000e-02, 1.000e-01, 6.000e-03, 1.500e-02, area=1.920e-03, Iyy=2.580e-06, Izz=1.820e-06, J=3.570e-06)
        S['120  x  40'] = {}
        S['120  x  40']['3.0'] = RHS(4.000e-02, 1.200e-01, 3.000e-03, 7.500e-03, area=9.010e-04, Iyy=1.480e-06, Izz=2.580e-07, J=7.460e-07)
        S['120  x  40']['4.0'] = RHS(4.000e-02, 1.200e-01, 4.000e-03, 1.000e-02, area=1.170e-03, Iyy=1.870e-06, Izz=3.190e-07, J=9.420e-07)
        S['120  x  40']['5.0'] = RHS(4.000e-02, 1.200e-01, 5.000e-03, 1.250e-02, area=1.440e-03, Iyy=2.210e-06, Izz=3.690e-07, J=1.110e-06)
        S['120  x  60'] = {}
        S['120  x  60']['3.0'] = RHS(6.000e-02, 1.200e-01, 3.000e-03, 7.500e-03, area=1.020e-03, Iyy=1.890e-06, Izz=6.440e-07, J=1.560e-06)
        S['120  x  60']['3.5'] = RHS(6.000e-02, 1.200e-01, 3.500e-03, 8.750e-03, area=1.180e-03, Iyy=2.160e-06, Izz=7.310e-07, J=1.790e-06)
        S['120  x  60']['4.0'] = RHS(6.000e-02, 1.200e-01, 4.000e-03, 1.000e-02, area=1.330e-03, Iyy=2.410e-06, Izz=8.120e-07, J=2.010e-06)
        S['120  x  60']['5.0'] = RHS(6.000e-02, 1.200e-01, 5.000e-03, 1.250e-02, area=1.640e-03, Iyy=2.870e-06, Izz=9.600e-07, J=2.420e-06)
        S['120  x  60']['6.0'] = RHS(6.000e-02, 1.200e-01, 6.000e-03, 1.500e-02, area=1.920e-03, Iyy=3.280e-06, Izz=1.090e-06, J=2.800e-06)
        S['120  x  80'] = {}
        S['120  x  80']['3.0'] = RHS(8.000e-02, 1.200e-01, 3.000e-03, 7.500e-03, area=1.140e-03, Iyy=2.300e-06, Izz=1.230e-06, J=2.550e-06)
        S['120  x  80']['4.0'] = RHS(8.000e-02, 1.200e-01, 4.000e-03, 1.000e-02, area=1.490e-03, Iyy=2.950e-06, Izz=1.570e-06, J=3.310e-06)
        S['120  x  80']['5.0'] = RHS(8.000e-02, 1.200e-01, 5.000e-03, 1.250e-02, area=1.840e-03, Iyy=3.530e-06, Izz=1.880e-06, J=4.020e-06)
        S['120  x  80']['6.0'] = RHS(8.000e-02, 1.200e-01, 6.000e-03, 1.500e-02, area=2.160e-03, Iyy=4.060e-06, Izz=2.150e-06, J=4.690e-06)
        S['120  x  80']['8.0'] = RHS(8.000e-02, 1.200e-01, 8.000e-03, 2.000e-02, area=2.720e-03, Iyy=4.760e-06, Izz=2.520e-06, J=5.840e-06)
        S['140  x  80'] = {}
        S['140  x  80']['3.0'] = RHS(8.000e-02, 1.400e-01, 3.000e-03, 7.500e-03, area=1.260e-03, Iyy=3.340e-06, Izz=1.410e-06, J=3.170e-06)
        S['140  x  80']['4.0'] = RHS(8.000e-02, 1.400e-01, 4.000e-03, 1.000e-02, area=1.650e-03, Iyy=4.300e-06, Izz=1.800e-06, J=4.120e-06)
        S['140  x  80']['5.0'] = RHS(8.000e-02, 1.400e-01, 5.000e-03, 1.250e-02, area=2.040e-03, Iyy=5.170e-06, Izz=2.160e-06, J=5.010e-06)
        S['140  x  80']['6.0'] = RHS(8.000e-02, 1.400e-01, 6.000e-03, 1.500e-02, area=2.400e-03, Iyy=5.970e-06, Izz=2.480e-06, J=5.840e-06)
        S['140  x  80']['8.0'] = RHS(8.000e-02, 1.400e-01, 8.000e-03, 2.000e-02, area=3.040e-03, Iyy=7.080e-06, Izz=2.930e-06, J=7.310e-06)
        S['140  x  80']['10.0'] = RHS(8.000e-02, 1.400e-01, 1.000e-02, 2.500e-02, area=3.660e-03, Iyy=8.040e-06, Izz=3.300e-06, J=8.510e-06)
        S['150  x  100'] = {}
        S['150  x  100']['3.0'] = RHS(1.000e-01, 1.500e-01, 3.000e-03, 7.500e-03, area=1.440e-03, Iyy=4.610e-06, Izz=2.480e-06, J=5.070e-06)
        S['150  x  100']['4.0'] = RHS(1.000e-01, 1.500e-01, 4.000e-03, 1.000e-02, area=1.890e-03, Iyy=5.950e-06, Izz=3.190e-06, J=6.620e-06)
        S['150  x  100']['5.0'] = RHS(1.000e-01, 1.500e-01, 5.000e-03, 1.250e-02, area=2.340e-03, Iyy=7.190e-06, Izz=3.840e-06, J=8.090e-06)
        S['150  x  100']['6.0'] = RHS(1.000e-01, 1.500e-01, 6.000e-03, 1.500e-02, area=2.760e-03, Iyy=8.350e-06, Izz=4.440e-06, J=9.480e-06)
        S['150  x  100']['8.0'] = RHS(1.000e-01, 1.500e-01, 8.000e-03, 2.000e-02, area=3.520e-03, Iyy=1.010e-05, Izz=5.360e-06, J=1.210e-05)
        S['150  x  100']['10.0'] = RHS(1.000e-01, 1.500e-01, 1.000e-02, 2.500e-02, area=4.260e-03, Iyy=1.160e-05, Izz=6.140e-06, J=1.430e-05)
        S['160  x  80'] = {}
        S['160  x  80']['3.0'] = RHS(8.000e-02, 1.600e-01, 3.000e-03, 7.500e-03, area=1.380e-03, Iyy=4.640e-06, Izz=1.590e-06, J=3.800e-06)
        S['160  x  80']['4.0'] = RHS(8.000e-02, 1.600e-01, 4.000e-03, 1.000e-02, area=1.810e-03, Iyy=5.980e-06, Izz=2.040e-06, J=4.940e-06)
        S['160  x  80']['5.0'] = RHS(8.000e-02, 1.600e-01, 5.000e-03, 1.250e-02, area=2.240e-03, Iyy=7.220e-06, Izz=2.440e-06, J=6.010e-06)
        S['160  x  80']['6.0'] = RHS(8.000e-02, 1.600e-01, 6.000e-03, 1.500e-02, area=2.640e-03, Iyy=8.360e-06, Izz=2.810e-06, J=7.020e-06)
        S['160  x  80']['8.0'] = RHS(8.000e-02, 1.600e-01, 8.000e-03, 2.000e-02, area=3.360e-03, Iyy=1.000e-05, Izz=3.350e-06, J=8.820e-06)
        S['160  x  80']['10.0'] = RHS(8.000e-02, 1.600e-01, 1.000e-02, 2.500e-02, area=4.060e-03, Iyy=1.150e-05, Izz=3.800e-06, J=1.030e-05)
        S['180  x  80'] = {}
        S['180  x  80']['3.0'] = RHS(8.000e-02, 1.800e-01, 3.000e-03, 7.500e-03, area=1.500e-03, Iyy=6.210e-06, Izz=1.770e-06, J=4.450e-06)
        S['180  x  80']['4.0'] = RHS(8.000e-02, 1.800e-01, 4.000e-03, 1.000e-02, area=1.970e-03, Iyy=8.020e-06, Izz=2.270e-06, J=5.780e-06)
        S['180  x  80']['5.0'] = RHS(8.000e-02, 1.800e-01, 5.000e-03, 1.250e-02, area=2.440e-03, Iyy=9.710e-06, Izz=2.720e-06, J=7.040e-06)
        S['180  x  80']['6.0'] = RHS(8.000e-02, 1.800e-01, 6.000e-03, 1.500e-02, area=2.880e-03, Iyy=1.130e-05, Izz=3.140e-06, J=8.230e-06)
        S['180  x  80']['8.0'] = RHS(8.000e-02, 1.800e-01, 8.000e-03, 2.000e-02, area=3.680e-03, Iyy=1.360e-05, Izz=3.770e-06, J=1.040e-05)
        S['180  x  80']['10.0'] = RHS(8.000e-02, 1.800e-01, 1.000e-02, 2.500e-02, area=4.460e-03, Iyy=1.570e-05, Izz=4.290e-06, J=1.210e-05)
        S['180  x  100'] = {}
        S['180  x  100']['4.0'] = RHS(1.000e-01, 1.800e-01, 4.000e-03, 1.000e-02, area=2.130e-03, Iyy=9.260e-06, Izz=3.740e-06, J=8.540e-06)
        S['180  x  100']['5.0'] = RHS(1.000e-01, 1.800e-01, 5.000e-03, 1.250e-02, area=2.640e-03, Iyy=1.120e-05, Izz=4.520e-06, J=1.040e-05)
        S['180  x  100']['6.0'] = RHS(1.000e-01, 1.800e-01, 6.000e-03, 1.500e-02, area=3.120e-03, Iyy=1.310e-05, Izz=5.240e-06, J=1.230e-05)
        S['180  x  100']['8.0'] = RHS(1.000e-01, 1.800e-01, 8.000e-03, 2.000e-02, area=4.000e-03, Iyy=1.600e-05, Izz=6.370e-06, J=1.560e-05)
        S['180  x  100']['10.0'] = RHS(1.000e-01, 1.800e-01, 1.000e-02, 2.500e-02, area=4.860e-03, Iyy=1.860e-05, Izz=7.360e-06, J=1.860e-05)
        S['200  x  100'] = {}
        S['200  x  100']['4.0'] = RHS(1.000e-01, 2.000e-01, 4.000e-03, 1.000e-02, area=2.290e-03, Iyy=1.200e-05, Izz=4.110e-06, J=9.850e-06)
        S['200  x  100']['5.0'] = RHS(1.000e-01, 2.000e-01, 5.000e-03, 1.250e-02, area=2.840e-03, Iyy=1.460e-05, Izz=4.970e-06, J=1.210e-05)
        S['200  x  100']['6.0'] = RHS(1.000e-01, 2.000e-01, 6.000e-03, 1.500e-02, area=3.360e-03, Iyy=1.700e-05, Izz=5.770e-06, J=1.420e-05)
        S['200  x  100']['8.0'] = RHS(1.000e-01, 2.000e-01, 8.000e-03, 2.000e-02, area=4.320e-03, Iyy=2.090e-05, Izz=7.050e-06, J=1.810e-05)
        S['200  x  100']['10.0'] = RHS(1.000e-01, 2.000e-01, 1.000e-02, 2.500e-02, area=5.260e-03, Iyy=2.440e-05, Izz=8.180e-06, J=2.150e-05)
        S['200  x  120'] = {}
        S['200  x  120']['4.0'] = RHS(1.200e-01, 2.000e-01, 4.000e-03, 1.000e-02, area=2.450e-03, Iyy=1.350e-05, Izz=6.180e-06, J=1.340e-05)
        S['200  x  120']['5.0'] = RHS(1.200e-01, 2.000e-01, 5.000e-03, 1.250e-02, area=3.040e-03, Iyy=1.650e-05, Izz=7.500e-06, J=1.650e-05)
        S['200  x  120']['6.0'] = RHS(1.200e-01, 2.000e-01, 6.000e-03, 1.500e-02, area=3.600e-03, Iyy=1.930e-05, Izz=8.740e-06, J=1.950e-05)
        S['200  x  120']['8.0'] = RHS(1.200e-01, 2.000e-01, 8.000e-03, 2.000e-02, area=4.640e-03, Iyy=2.390e-05, Izz=1.080e-05, J=2.510e-05)
        S['200  x  120']['10.0'] = RHS(1.200e-01, 2.000e-01, 1.000e-02, 2.500e-02, area=5.660e-03, Iyy=2.810e-05, Izz=1.260e-05, J=3.010e-05)
        S['200  x  150'] = {}
        S['200  x  150']['4.0'] = RHS(1.500e-01, 2.000e-01, 4.000e-03, 1.000e-02, area=2.690e-03, Iyy=1.580e-05, Izz=1.020e-05, J=1.940e-05)
        S['200  x  150']['5.0'] = RHS(1.500e-01, 2.000e-01, 5.000e-03, 1.250e-02, area=3.340e-03, Iyy=1.940e-05, Izz=1.240e-05, J=2.390e-05)
        S['200  x  150']['6.0'] = RHS(1.500e-01, 2.000e-01, 6.000e-03, 1.500e-02, area=3.960e-03, Iyy=2.270e-05, Izz=1.460e-05, J=2.830e-05)
        S['200  x  150']['8.0'] = RHS(1.500e-01, 2.000e-01, 8.000e-03, 2.000e-02, area=5.120e-03, Iyy=2.830e-05, Izz=1.820e-05, J=3.660e-05)
        S['200  x  150']['10.0'] = RHS(1.500e-01, 2.000e-01, 1.000e-02, 2.500e-02, area=6.260e-03, Iyy=3.350e-05, Izz=2.140e-05, J=4.430e-05)
        S['250  x  150'] = {}
        S['250  x  150']['5.0'] = RHS(1.500e-01, 2.500e-01, 5.000e-03, 1.250e-02, area=3.840e-03, Iyy=3.300e-05, Izz=1.510e-05, J=3.280e-05)
        S['250  x  150']['6.0'] = RHS(1.500e-01, 2.500e-01, 6.000e-03, 1.500e-02, area=4.560e-03, Iyy=3.890e-05, Izz=1.770e-05, J=3.890e-05)
        S['250  x  150']['6.3'] = RHS(1.500e-01, 2.500e-01, 6.300e-03, 1.575e-02, area=4.740e-03, Iyy=4.000e-05, Izz=1.820e-05, J=4.080e-05)
        S['250  x  150']['8.0'] = RHS(1.500e-01, 2.500e-01, 8.000e-03, 2.000e-02, area=5.920e-03, Iyy=4.890e-05, Izz=2.220e-05, J=5.050e-05)
        S['250  x  150']['10.0'] = RHS(1.500e-01, 2.500e-01, 1.000e-02, 2.500e-02, area=7.260e-03, Iyy=5.820e-05, Izz=2.630e-05, J=6.120e-05)
        S['250  x  150']['12.0'] = RHS(1.500e-01, 2.500e-01, 1.200e-02, 3.000e-02, area=8.410e-03, Iyy=6.460e-05, Izz=2.920e-05, J=7.090e-05)
        S['250  x  150']['12.5'] = RHS(1.500e-01, 2.500e-01, 1.250e-02, 3.125e-02, area=8.700e-03, Iyy=6.630e-05, Izz=3.000e-05, J=7.320e-05)
        S['300  x  100'] = {}
        S['300  x  100']['6.0'] = RHS(1.000e-01, 3.000e-01, 6.000e-03, 1.500e-02, area=4.560e-03, Iyy=4.780e-05, Izz=8.420e-06, J=2.400e-05)
        S['300  x  100']['8.0'] = RHS(1.000e-01, 3.000e-01, 8.000e-03, 2.000e-02, area=5.920e-03, Iyy=5.980e-05, Izz=1.040e-05, J=3.080e-05)
        S['300  x  100']['10.0'] = RHS(1.000e-01, 3.000e-01, 1.000e-02, 2.500e-02, area=7.260e-03, Iyy=7.110e-05, Izz=1.220e-05, J=3.680e-05)
        S['300  x  200'] = {}
        S['300  x  200']['6.0'] = RHS(2.000e-01, 3.000e-01, 6.000e-03, 1.500e-02, area=5.760e-03, Iyy=7.370e-05, Izz=3.960e-05, J=8.120e-05)
        S['300  x  200']['6.3'] = RHS(2.000e-01, 3.000e-01, 6.300e-03, 1.575e-02, area=6.000e-03, Iyy=7.620e-05, Izz=4.100e-05, J=8.520e-05)
        S['300  x  200']['8.0'] = RHS(2.000e-01, 3.000e-01, 8.000e-03, 2.000e-02, area=7.520e-03, Iyy=9.390e-05, Izz=5.040e-05, J=1.060e-04)
        S['300  x  200']['10.0'] = RHS(2.000e-01, 3.000e-01, 1.000e-02, 2.500e-02, area=9.260e-03, Iyy=1.130e-04, Izz=6.060e-05, J=1.300e-04)
        S['300  x  200']['12.0'] = RHS(2.000e-01, 3.000e-01, 1.200e-02, 3.000e-02, area=1.080e-02, Iyy=1.280e-04, Izz=6.850e-05, J=1.520e-04)
        S['300  x  200']['12.5'] = RHS(2.000e-01, 3.000e-01, 1.250e-02, 3.125e-02, area=1.120e-02, Iyy=1.320e-04, Izz=7.060e-05, J=1.580e-04)
        S['400  x  200'] = {}
        S['400  x  200']['6.0'] = RHS(2.000e-01, 4.000e-01, 6.000e-03, 1.500e-02, area=6.960e-03, Iyy=1.480e-04, Izz=5.090e-05, J=1.210e-04)
        S['400  x  200']['6.3'] = RHS(2.000e-01, 4.000e-01, 6.300e-03, 1.575e-02, area=7.260e-03, Iyy=1.530e-04, Izz=5.290e-05, J=1.270e-04)
        S['400  x  200']['8.0'] = RHS(2.000e-01, 4.000e-01, 8.000e-03, 2.000e-02, area=9.120e-03, Iyy=1.900e-04, Izz=6.520e-05, J=1.580e-04)
        S['400  x  200']['10.0'] = RHS(2.000e-01, 4.000e-01, 1.000e-02, 2.500e-02, area=1.130e-02, Iyy=2.300e-04, Izz=7.860e-05, J=1.940e-04)
        S['400  x  200']['12.0'] = RHS(2.000e-01, 4.000e-01, 1.200e-02, 3.000e-02, area=1.320e-02, Iyy=2.620e-04, Izz=8.980e-05, J=2.280e-04)
        S['400  x  200']['12.5'] = RHS(2.000e-01, 4.000e-01, 1.250e-02, 3.125e-02, area=1.370e-02, Iyy=2.710e-04, Izz=9.260e-05, J=2.360e-04)
        S['450  x  250'] = {}
        S['450  x  250']['6.0'] = RHS(2.500e-01, 4.500e-01, 6.000e-03, 1.500e-02, area=8.160e-03, Iyy=2.270e-04, Izz=9.240e-05, J=2.070e-04)
        S['450  x  250']['6.3'] = RHS(2.500e-01, 4.500e-01, 6.300e-03, 1.575e-02, area=8.520e-03, Iyy=2.360e-04, Izz=9.620e-05, J=2.170e-04)
        S['450  x  250']['8.0'] = RHS(2.500e-01, 4.500e-01, 8.000e-03, 2.000e-02, area=1.070e-02, Iyy=2.930e-04, Izz=1.190e-04, J=2.720e-04)
        S['450  x  250']['10.0'] = RHS(2.500e-01, 4.500e-01, 1.000e-02, 2.500e-02, area=1.330e-02, Iyy=3.570e-04, Izz=1.450e-04, J=3.350e-04)
        S['450  x  250']['12.0'] = RHS(2.500e-01, 4.500e-01, 1.200e-02, 3.000e-02, area=1.560e-02, Iyy=4.110e-04, Izz=1.670e-04, J=3.960e-04)
        S['450  x  250']['12.5'] = RHS(2.500e-01, 4.500e-01, 1.250e-02, 3.125e-02, area=1.620e-02, Iyy=4.250e-04, Izz=1.720e-04, J=4.110e-04)
        S['500  x  300'] = {}
        S['500  x  300']['6.0'] = RHS(3.000e-01, 5.000e-01, 6.000e-03, 1.500e-02, area=9.360e-03, Iyy=3.300e-04, Izz=1.520e-04, J=3.240e-04)
        S['500  x  300']['6.3'] = RHS(3.000e-01, 5.000e-01, 6.300e-03, 1.575e-02, area=9.780e-03, Iyy=3.430e-04, Izz=1.580e-04, J=3.410e-04)
        S['500  x  300']['8.0'] = RHS(3.000e-01, 5.000e-01, 8.000e-03, 2.000e-02, area=1.230e-02, Iyy=4.280e-04, Izz=1.960e-04, J=4.280e-04)
        S['500  x  300']['10.0'] = RHS(3.000e-01, 5.000e-01, 1.000e-02, 2.500e-02, area=1.530e-02, Iyy=5.230e-04, Izz=2.390e-04, J=5.270e-04)
        S['500  x  300']['12.0'] = RHS(3.000e-01, 5.000e-01, 1.200e-02, 3.000e-02, area=1.800e-02, Iyy=6.060e-04, Izz=2.770e-04, J=6.260e-04)
        S['500  x  300']['12.5'] = RHS(3.000e-01, 5.000e-01, 1.250e-02, 3.125e-02, area=1.870e-02, Iyy=6.270e-04, Izz=2.870e-04, J=6.500e-04)

        return S # No. data rows = 137

    @staticmethod
    def RHS_ColdFormed():
        if BlueBook.__RHS_ColdFormed is None:
            BlueBook.__RHS_ColdFormed = BlueBook.__create_RHS_ColdFormed()
        return BlueBook.__RHS_ColdFormed

    __RHS_HotFinished = None

    @staticmethod
    def __create_RHS_HotFinished():
        S = {}
        S['50  x  30'] = {}
        S['50  x  30']['2.9'] = RHS(3.000e-02, 5.000e-02, 2.900e-03, 4.350e-03, area=4.210e-04, Iyy=1.320e-07, Izz=5.800e-08, J=1.320e-07)
        S['50  x  30']['3.0'] = RHS(3.000e-02, 5.000e-02, 3.000e-03, 4.500e-03, area=4.340e-04, Iyy=1.360e-07, Izz=5.940e-08, J=1.350e-07)
        S['50  x  30']['3.2'] = RHS(3.000e-02, 5.000e-02, 3.200e-03, 4.800e-03, area=4.600e-04, Iyy=1.420e-07, Izz=6.200e-08, J=1.420e-07)
        S['50  x  30']['3.6'] = RHS(3.000e-02, 5.000e-02, 3.600e-03, 5.400e-03, area=5.100e-04, Iyy=1.540e-07, Izz=6.670e-08, J=1.540e-07)
        S['50  x  30']['4.0'] = RHS(3.000e-02, 5.000e-02, 4.000e-03, 6.000e-03, area=5.590e-04, Iyy=1.650e-07, Izz=7.080e-08, J=1.660e-07)
        S['50  x  30']['5.0'] = RHS(3.000e-02, 5.000e-02, 5.000e-03, 7.500e-03, area=6.730e-04, Iyy=1.870e-07, Izz=7.890e-08, J=1.900e-07)
        S['50  x  30']['5.6'] = RHS(3.000e-02, 5.000e-02, 5.600e-03, 8.400e-03, area=7.370e-04, Iyy=1.970e-07, Izz=8.230e-08, J=2.010e-07)
        S['50  x  30']['6.3'] = RHS(3.000e-02, 5.000e-02, 6.300e-03, 9.450e-03, area=8.070e-04, Iyy=2.060e-07, Izz=8.500e-08, J=2.110e-07)
        S['60  x  40'] = {}
        S['60  x  40']['2.9'] = RHS(4.000e-02, 6.000e-02, 2.900e-03, 4.350e-03, area=5.370e-04, Iyy=2.580e-07, Izz=1.350e-07, J=2.840e-07)
        S['60  x  40']['3.0'] = RHS(4.000e-02, 6.000e-02, 3.000e-03, 4.500e-03, area=5.540e-04, Iyy=2.650e-07, Izz=1.390e-07, J=2.920e-07)
        S['60  x  40']['3.2'] = RHS(4.000e-02, 6.000e-02, 3.200e-03, 4.800e-03, area=5.880e-04, Iyy=2.780e-07, Izz=1.460e-07, J=3.080e-07)
        S['60  x  40']['3.6'] = RHS(4.000e-02, 6.000e-02, 3.600e-03, 5.400e-03, area=6.540e-04, Iyy=3.040e-07, Izz=1.590e-07, J=3.380e-07)
        S['60  x  40']['4.0'] = RHS(4.000e-02, 6.000e-02, 4.000e-03, 6.000e-03, area=7.190e-04, Iyy=3.280e-07, Izz=1.700e-07, J=3.670e-07)
        S['60  x  40']['5.0'] = RHS(4.000e-02, 6.000e-02, 5.000e-03, 7.500e-03, area=8.730e-04, Iyy=3.810e-07, Izz=1.950e-07, J=4.300e-07)
        S['60  x  40']['5.6'] = RHS(4.000e-02, 6.000e-02, 5.600e-03, 8.400e-03, area=9.610e-04, Iyy=4.070e-07, Izz=2.070e-07, J=4.620e-07)
        S['60  x  40']['6.3'] = RHS(4.000e-02, 6.000e-02, 6.300e-03, 9.450e-03, area=1.060e-03, Iyy=4.340e-07, Izz=2.190e-07, J=4.950e-07)
        S['60  x  40']['7.1'] = RHS(4.000e-02, 6.000e-02, 7.100e-03, 1.065e-02, area=1.160e-03, Iyy=4.590e-07, Izz=2.290e-07, J=5.270e-07)
        S['60  x  40']['8.0'] = RHS(4.000e-02, 6.000e-02, 8.000e-03, 1.200e-02, area=1.280e-03, Iyy=4.790e-07, Izz=2.370e-07, J=5.540e-07)
        S['80  x  40'] = {}
        S['80  x  40']['2.9'] = RHS(4.000e-02, 8.000e-02, 2.900e-03, 4.350e-03, area=6.530e-04, Iyy=5.270e-07, Izz=1.750e-07, J=4.260e-07)
        S['80  x  40']['3.0'] = RHS(4.000e-02, 8.000e-02, 3.000e-03, 4.500e-03, area=6.740e-04, Iyy=5.420e-07, Izz=1.800e-07, J=4.380e-07)
        S['80  x  40']['3.2'] = RHS(4.000e-02, 8.000e-02, 3.200e-03, 4.800e-03, area=7.160e-04, Iyy=5.720e-07, Izz=1.890e-07, J=4.620e-07)
        S['80  x  40']['3.6'] = RHS(4.000e-02, 8.000e-02, 3.600e-03, 5.400e-03, area=7.980e-04, Iyy=6.280e-07, Izz=2.060e-07, J=5.080e-07)
        S['80  x  40']['4.0'] = RHS(4.000e-02, 8.000e-02, 4.000e-03, 6.000e-03, area=8.790e-04, Iyy=6.820e-07, Izz=2.220e-07, J=5.520e-07)
        S['80  x  40']['5.0'] = RHS(4.000e-02, 8.000e-02, 5.000e-03, 7.500e-03, area=1.070e-03, Iyy=8.030e-07, Izz=2.570e-07, J=6.510e-07)
        S['80  x  40']['5.6'] = RHS(4.000e-02, 8.000e-02, 5.600e-03, 8.400e-03, area=1.180e-03, Iyy=8.670e-07, Izz=2.740e-07, J=7.020e-07)
        S['80  x  40']['6.3'] = RHS(4.000e-02, 8.000e-02, 6.300e-03, 9.450e-03, area=1.310e-03, Iyy=9.330e-07, Izz=2.920e-07, J=7.560e-07)
        S['80  x  40']['7.1'] = RHS(4.000e-02, 8.000e-02, 7.100e-03, 1.065e-02, area=1.450e-03, Iyy=9.980e-07, Izz=3.070e-07, J=8.090e-07)
        S['80  x  40']['8.0'] = RHS(4.000e-02, 8.000e-02, 8.000e-03, 1.200e-02, area=1.600e-03, Iyy=1.060e-06, Izz=3.210e-07, J=8.580e-07)
        S['90  x  50'] = {}
        S['90  x  50']['3.0'] = RHS(5.000e-02, 9.000e-02, 3.000e-03, 4.500e-03, area=7.940e-04, Iyy=8.440e-07, Izz=3.350e-07, J=7.650e-07)
        S['90  x  50']['3.2'] = RHS(5.000e-02, 9.000e-02, 3.200e-03, 4.800e-03, area=8.440e-04, Iyy=8.910e-07, Izz=3.530e-07, J=8.090e-07)
        S['90  x  50']['3.6'] = RHS(5.000e-02, 9.000e-02, 3.600e-03, 5.400e-03, area=9.420e-04, Iyy=9.830e-07, Izz=3.870e-07, J=8.940e-07)
        S['90  x  50']['4.0'] = RHS(5.000e-02, 9.000e-02, 4.000e-03, 6.000e-03, area=1.040e-03, Iyy=1.070e-06, Izz=4.190e-07, J=9.750e-07)
        S['90  x  50']['5.0'] = RHS(5.000e-02, 9.000e-02, 5.000e-03, 7.500e-03, area=1.270e-03, Iyy=1.270e-06, Izz=4.920e-07, J=1.160e-06)
        S['90  x  50']['5.6'] = RHS(5.000e-02, 9.000e-02, 5.600e-03, 8.400e-03, area=1.410e-03, Iyy=1.380e-06, Izz=5.300e-07, J=1.270e-06)
        S['90  x  50']['6.3'] = RHS(5.000e-02, 9.000e-02, 6.300e-03, 9.450e-03, area=1.560e-03, Iyy=1.500e-06, Izz=5.700e-07, J=1.380e-06)
        S['90  x  50']['7.1'] = RHS(5.000e-02, 9.000e-02, 7.100e-03, 1.065e-02, area=1.730e-03, Iyy=1.620e-06, Izz=6.090e-07, J=1.490e-06)
        S['90  x  50']['8.0'] = RHS(5.000e-02, 9.000e-02, 8.000e-03, 1.200e-02, area=1.920e-03, Iyy=1.740e-06, Izz=6.460e-07, J=1.600e-06)
        S['90  x  50']['8.8'] = RHS(5.000e-02, 9.000e-02, 8.800e-03, 1.320e-02, area=2.070e-03, Iyy=1.830e-06, Izz=6.720e-07, J=1.690e-06)
        S['100  x  50'] = {}
        S['100  x  50']['3.2'] = RHS(5.000e-02, 1.000e-01, 3.200e-03, 4.800e-03, area=9.080e-04, Iyy=1.160e-06, Izz=3.880e-07, J=9.340e-07)
        S['100  x  50']['3.6'] = RHS(5.000e-02, 1.000e-01, 3.600e-03, 5.400e-03, area=1.010e-03, Iyy=1.280e-06, Izz=4.260e-07, J=1.030e-06)
        S['100  x  50']['4.0'] = RHS(5.000e-02, 1.000e-01, 4.000e-03, 6.000e-03, area=1.120e-03, Iyy=1.400e-06, Izz=4.620e-07, J=1.130e-06)
        S['100  x  50']['5.0'] = RHS(5.000e-02, 1.000e-01, 5.000e-03, 7.500e-03, area=1.370e-03, Iyy=1.670e-06, Izz=5.430e-07, J=1.350e-06)
        S['100  x  50']['5.6'] = RHS(5.000e-02, 1.000e-01, 5.600e-03, 8.400e-03, area=1.520e-03, Iyy=1.810e-06, Izz=5.860e-07, J=1.470e-06)
        S['100  x  50']['6.3'] = RHS(5.000e-02, 1.000e-01, 6.300e-03, 9.450e-03, area=1.690e-03, Iyy=1.970e-06, Izz=6.300e-07, J=1.600e-06)
        S['100  x  50']['7.1'] = RHS(5.000e-02, 1.000e-01, 7.100e-03, 1.065e-02, area=1.870e-03, Iyy=2.140e-06, Izz=6.750e-07, J=1.730e-06)
        S['100  x  50']['8.0'] = RHS(5.000e-02, 1.000e-01, 8.000e-03, 1.200e-02, area=2.080e-03, Iyy=2.300e-06, Izz=7.170e-07, J=1.860e-06)
        S['100  x  50']['8.8'] = RHS(5.000e-02, 1.000e-01, 8.800e-03, 1.320e-02, area=2.250e-03, Iyy=2.430e-06, Izz=7.480e-07, J=1.970e-06)
        S['100  x  50']['10.0'] = RHS(5.000e-02, 1.000e-01, 1.000e-02, 1.500e-02, area=2.490e-03, Iyy=2.590e-06, Izz=7.840e-07, J=2.090e-06)
        S['100  x  60'] = {}
        S['100  x  60']['3.2'] = RHS(6.000e-02, 1.000e-01, 3.200e-03, 4.800e-03, area=9.720e-04, Iyy=1.310e-06, Izz=5.880e-07, J=1.290e-06)
        S['100  x  60']['3.6'] = RHS(6.000e-02, 1.000e-01, 3.600e-03, 5.400e-03, area=1.090e-03, Iyy=1.450e-06, Izz=6.480e-07, J=1.420e-06)
        S['100  x  60']['4.0'] = RHS(6.000e-02, 1.000e-01, 4.000e-03, 6.000e-03, area=1.200e-03, Iyy=1.580e-06, Izz=7.050e-07, J=1.560e-06)
        S['100  x  60']['5.0'] = RHS(6.000e-02, 1.000e-01, 5.000e-03, 7.500e-03, area=1.470e-03, Iyy=1.890e-06, Izz=8.360e-07, J=1.880e-06)
        S['100  x  60']['5.6'] = RHS(6.000e-02, 1.000e-01, 5.600e-03, 8.400e-03, area=1.630e-03, Iyy=2.060e-06, Izz=9.060e-07, J=2.050e-06)
        S['100  x  60']['6.3'] = RHS(6.000e-02, 1.000e-01, 6.300e-03, 9.450e-03, area=1.810e-03, Iyy=2.250e-06, Izz=9.810e-07, J=2.240e-06)
        S['100  x  60']['7.1'] = RHS(6.000e-02, 1.000e-01, 7.100e-03, 1.065e-02, area=2.020e-03, Iyy=2.440e-06, Izz=1.060e-06, J=2.450e-06)
        S['100  x  60']['8.0'] = RHS(6.000e-02, 1.000e-01, 8.000e-03, 1.200e-02, area=2.240e-03, Iyy=2.640e-06, Izz=1.130e-06, J=2.650e-06)
        S['100  x  60']['8.8'] = RHS(6.000e-02, 1.000e-01, 8.800e-03, 1.320e-02, area=2.420e-03, Iyy=2.790e-06, Izz=1.190e-06, J=2.820e-06)
        S['100  x  60']['10.0'] = RHS(6.000e-02, 1.000e-01, 1.000e-02, 1.500e-02, area=2.690e-03, Iyy=2.990e-06, Izz=1.260e-06, J=3.040e-06)
        S['120  x  60'] = {}
        S['120  x  60']['3.6'] = RHS(6.000e-02, 1.200e-01, 3.600e-03, 5.400e-03, area=1.230e-03, Iyy=2.270e-06, Izz=7.630e-07, J=1.830e-06)
        S['120  x  60']['4.0'] = RHS(6.000e-02, 1.200e-01, 4.000e-03, 6.000e-03, area=1.360e-03, Iyy=2.490e-06, Izz=8.310e-07, J=2.010e-06)
        S['120  x  60']['5.0'] = RHS(6.000e-02, 1.200e-01, 5.000e-03, 7.500e-03, area=1.670e-03, Iyy=2.990e-06, Izz=9.880e-07, J=2.420e-06)
        S['120  x  60']['5.6'] = RHS(6.000e-02, 1.200e-01, 5.600e-03, 8.400e-03, area=1.860e-03, Iyy=3.270e-06, Izz=1.070e-06, J=2.650e-06)
        S['120  x  60']['6.3'] = RHS(6.000e-02, 1.200e-01, 6.300e-03, 9.450e-03, area=2.070e-03, Iyy=3.580e-06, Izz=1.160e-06, J=2.900e-06)
        S['120  x  60']['7.1'] = RHS(6.000e-02, 1.200e-01, 7.100e-03, 1.065e-02, area=2.300e-03, Iyy=3.910e-06, Izz=1.260e-06, J=3.170e-06)
        S['120  x  60']['8.0'] = RHS(6.000e-02, 1.200e-01, 8.000e-03, 1.200e-02, area=2.560e-03, Iyy=4.250e-06, Izz=1.350e-06, J=3.440e-06)
        S['120  x  60']['8.8'] = RHS(6.000e-02, 1.200e-01, 8.800e-03, 1.320e-02, area=2.780e-03, Iyy=4.520e-06, Izz=1.420e-06, J=3.660e-06)
        S['120  x  60']['10.0'] = RHS(6.000e-02, 1.200e-01, 1.000e-02, 1.500e-02, area=3.090e-03, Iyy=4.880e-06, Izz=1.520e-06, J=3.960e-06)
        S['120  x  80'] = {}
        S['120  x  80']['3.6'] = RHS(8.000e-02, 1.200e-01, 3.600e-03, 5.400e-03, area=1.370e-03, Iyy=2.760e-06, Izz=1.470e-06, J=3.010e-06)
        S['120  x  80']['4.0'] = RHS(8.000e-02, 1.200e-01, 4.000e-03, 6.000e-03, area=1.520e-03, Iyy=3.030e-06, Izz=1.610e-06, J=3.300e-06)
        S['120  x  80']['5.0'] = RHS(8.000e-02, 1.200e-01, 5.000e-03, 7.500e-03, area=1.870e-03, Iyy=3.650e-06, Izz=1.930e-06, J=4.010e-06)
        S['120  x  80']['5.6'] = RHS(8.000e-02, 1.200e-01, 5.600e-03, 8.400e-03, area=2.080e-03, Iyy=4.010e-06, Izz=2.110e-06, J=4.420e-06)
        S['120  x  80']['6.3'] = RHS(8.000e-02, 1.200e-01, 6.300e-03, 9.450e-03, area=2.320e-03, Iyy=4.400e-06, Izz=2.300e-06, J=4.870e-06)
        S['120  x  80']['7.1'] = RHS(8.000e-02, 1.200e-01, 7.100e-03, 1.065e-02, area=2.580e-03, Iyy=4.820e-06, Izz=2.510e-06, J=5.350e-06)
        S['120  x  80']['8.0'] = RHS(8.000e-02, 1.200e-01, 8.000e-03, 1.200e-02, area=2.880e-03, Iyy=5.250e-06, Izz=2.730e-06, J=5.870e-06)
        S['120  x  80']['8.8'] = RHS(8.000e-02, 1.200e-01, 8.800e-03, 1.320e-02, area=3.130e-03, Iyy=5.610e-06, Izz=2.900e-06, J=6.290e-06)
        S['120  x  80']['10.0'] = RHS(8.000e-02, 1.200e-01, 1.000e-02, 1.500e-02, area=3.490e-03, Iyy=6.090e-06, Izz=3.130e-06, J=6.880e-06)
        S['150  x  100'] = {}
        S['150  x  100']['4.0'] = RHS(1.000e-01, 1.500e-01, 4.000e-03, 6.000e-03, area=1.920e-03, Iyy=6.070e-06, Izz=3.240e-06, J=6.600e-06)
        S['150  x  100']['5.0'] = RHS(1.000e-01, 1.500e-01, 5.000e-03, 7.500e-03, area=2.370e-03, Iyy=7.390e-06, Izz=3.920e-06, J=8.070e-06)
        S['150  x  100']['5.6'] = RHS(1.000e-01, 1.500e-01, 5.600e-03, 8.400e-03, area=2.640e-03, Iyy=8.140e-06, Izz=4.310e-06, J=8.910e-06)
        S['150  x  100']['6.3'] = RHS(1.000e-01, 1.500e-01, 6.300e-03, 9.450e-03, area=2.950e-03, Iyy=8.980e-06, Izz=4.740e-06, J=9.860e-06)
        S['150  x  100']['7.1'] = RHS(1.000e-01, 1.500e-01, 7.100e-03, 1.065e-02, area=3.290e-03, Iyy=9.900e-06, Izz=5.200e-06, J=1.090e-05)
        S['150  x  100']['8.0'] = RHS(1.000e-01, 1.500e-01, 8.000e-03, 1.200e-02, area=3.680e-03, Iyy=1.090e-05, Izz=5.690e-06, J=1.200e-05)
        S['150  x  100']['8.8'] = RHS(1.000e-01, 1.500e-01, 8.800e-03, 1.320e-02, area=4.010e-03, Iyy=1.170e-05, Izz=6.100e-06, J=1.300e-05)
        S['150  x  100']['10.0'] = RHS(1.000e-01, 1.500e-01, 1.000e-02, 1.500e-02, area=4.490e-03, Iyy=1.280e-05, Izz=6.650e-06, J=1.430e-05)
        S['150  x  100']['11.0'] = RHS(1.000e-01, 1.500e-01, 1.100e-02, 1.650e-02, area=4.890e-03, Iyy=1.370e-05, Izz=7.070e-06, J=1.540e-05)
        S['150  x  100']['12.5'] = RHS(1.000e-01, 1.500e-01, 1.250e-02, 1.875e-02, area=5.460e-03, Iyy=1.490e-05, Izz=7.630e-06, J=1.680e-05)
        S['160  x  80'] = {}
        S['160  x  80']['4.0'] = RHS(8.000e-02, 1.600e-01, 4.000e-03, 6.000e-03, area=1.840e-03, Iyy=6.120e-06, Izz=2.070e-06, J=4.930e-06)
        S['160  x  80']['5.0'] = RHS(8.000e-02, 1.600e-01, 5.000e-03, 7.500e-03, area=2.270e-03, Iyy=7.440e-06, Izz=2.490e-06, J=6.000e-06)
        S['160  x  80']['5.6'] = RHS(8.000e-02, 1.600e-01, 5.600e-03, 8.400e-03, area=2.530e-03, Iyy=8.190e-06, Izz=2.730e-06, J=6.610e-06)
        S['160  x  80']['6.3'] = RHS(8.000e-02, 1.600e-01, 6.300e-03, 9.450e-03, area=2.820e-03, Iyy=9.030e-06, Izz=2.990e-06, J=7.300e-06)
        S['160  x  80']['7.1'] = RHS(8.000e-02, 1.600e-01, 7.100e-03, 1.065e-02, area=3.150e-03, Iyy=9.940e-06, Izz=3.270e-06, J=8.040e-06)
        S['160  x  80']['8.0'] = RHS(8.000e-02, 1.600e-01, 8.000e-03, 1.200e-02, area=3.520e-03, Iyy=1.090e-05, Izz=3.560e-06, J=8.830e-06)
        S['160  x  80']['8.8'] = RHS(8.000e-02, 1.600e-01, 8.800e-03, 1.320e-02, area=3.830e-03, Iyy=1.170e-05, Izz=3.790e-06, J=9.490e-06)
        S['160  x  80']['10.0'] = RHS(8.000e-02, 1.600e-01, 1.000e-02, 1.500e-02, area=4.290e-03, Iyy=1.280e-05, Izz=4.110e-06, J=1.040e-05)
        S['160  x  80']['11.0'] = RHS(8.000e-02, 1.600e-01, 1.100e-02, 1.650e-02, area=4.670e-03, Iyy=1.370e-05, Izz=4.350e-06, J=1.110e-05)
        S['160  x  80']['12.5'] = RHS(8.000e-02, 1.600e-01, 1.250e-02, 1.875e-02, area=5.210e-03, Iyy=1.480e-05, Izz=4.650e-06, J=1.200e-05)
        S['180  x  60'] = {}
        S['180  x  60']['4.0'] = RHS(6.000e-02, 1.800e-01, 4.000e-03, 6.000e-03, area=1.840e-03, Iyy=6.970e-06, Izz=1.210e-06, J=3.410e-06)
        S['180  x  60']['5.0'] = RHS(6.000e-02, 1.800e-01, 5.000e-03, 7.500e-03, area=2.270e-03, Iyy=8.460e-06, Izz=1.440e-06, J=4.110e-06)
        S['180  x  60']['5.6'] = RHS(6.000e-02, 1.800e-01, 5.600e-03, 8.400e-03, area=2.530e-03, Iyy=9.320e-06, Izz=1.570e-06, J=4.510e-06)
        S['180  x  60']['6.3'] = RHS(6.000e-02, 1.800e-01, 6.300e-03, 9.450e-03, area=2.820e-03, Iyy=1.030e-05, Izz=1.710e-06, J=4.950e-06)
        S['180  x  60']['7.1'] = RHS(6.000e-02, 1.800e-01, 7.100e-03, 1.065e-02, area=3.150e-03, Iyy=1.130e-05, Izz=1.860e-06, J=5.420e-06)
        S['180  x  60']['8.0'] = RHS(6.000e-02, 1.800e-01, 8.000e-03, 1.200e-02, area=3.520e-03, Iyy=1.240e-05, Izz=2.010e-06, J=5.900e-06)
        S['180  x  60']['8.8'] = RHS(6.000e-02, 1.800e-01, 8.800e-03, 1.320e-02, area=3.830e-03, Iyy=1.330e-05, Izz=2.120e-06, J=6.300e-06)
        S['180  x  60']['10.0'] = RHS(6.000e-02, 1.800e-01, 1.000e-02, 1.500e-02, area=4.290e-03, Iyy=1.460e-05, Izz=2.280e-06, J=6.830e-06)
        S['180  x  60']['11.0'] = RHS(6.000e-02, 1.800e-01, 1.100e-02, 1.650e-02, area=4.670e-03, Iyy=1.550e-05, Izz=2.380e-06, J=7.220e-06)
        S['180  x  60']['12.5'] = RHS(6.000e-02, 1.800e-01, 1.250e-02, 1.875e-02, area=5.210e-03, Iyy=1.680e-05, Izz=2.510e-06, J=7.700e-06)
        S['180  x  100'] = {}
        S['180  x  100']['4.0'] = RHS(1.000e-01, 1.800e-01, 4.000e-03, 6.000e-03, area=2.160e-03, Iyy=9.450e-06, Izz=3.790e-06, J=8.520e-06)
        S['180  x  100']['5.0'] = RHS(1.000e-01, 1.800e-01, 5.000e-03, 7.500e-03, area=2.670e-03, Iyy=1.150e-05, Izz=4.600e-06, J=1.040e-05)
        S['180  x  100']['5.6'] = RHS(1.000e-01, 1.800e-01, 5.600e-03, 8.400e-03, area=2.980e-03, Iyy=1.270e-05, Izz=5.060e-06, J=1.150e-05)
        S['180  x  100']['6.3'] = RHS(1.000e-01, 1.800e-01, 6.300e-03, 9.450e-03, area=3.330e-03, Iyy=1.410e-05, Izz=5.570e-06, J=1.280e-05)
        S['180  x  100']['7.1'] = RHS(1.000e-01, 1.800e-01, 7.100e-03, 1.065e-02, area=3.720e-03, Iyy=1.560e-05, Izz=6.130e-06, J=1.410e-05)
        S['180  x  100']['8.0'] = RHS(1.000e-01, 1.800e-01, 8.000e-03, 1.200e-02, area=4.160e-03, Iyy=1.710e-05, Izz=6.710e-06, J=1.560e-05)
        S['180  x  100']['8.8'] = RHS(1.000e-01, 1.800e-01, 8.800e-03, 1.320e-02, area=4.540e-03, Iyy=1.850e-05, Izz=7.200e-06, J=1.680e-05)
        S['180  x  100']['10.0'] = RHS(1.000e-01, 1.800e-01, 1.000e-02, 1.500e-02, area=5.090e-03, Iyy=2.040e-05, Izz=7.870e-06, J=1.860e-05)
        S['180  x  100']['11.0'] = RHS(1.000e-01, 1.800e-01, 1.100e-02, 1.650e-02, area=5.550e-03, Iyy=2.180e-05, Izz=8.390e-06, J=2.000e-05)
        S['180  x  100']['12.5'] = RHS(1.000e-01, 1.800e-01, 1.250e-02, 1.875e-02, area=6.210e-03, Iyy=2.380e-05, Izz=9.080e-06, J=2.190e-05)
        S['180  x  100']['14.2'] = RHS(1.000e-01, 1.800e-01, 1.420e-02, 2.130e-02, area=6.930e-03, Iyy=2.590e-05, Izz=9.740e-06, J=2.380e-05)
        S['200  x  100'] = {}
        S['200  x  100']['4.0'] = RHS(1.000e-01, 2.000e-01, 4.000e-03, 6.000e-03, area=2.320e-03, Iyy=1.220e-05, Izz=4.160e-06, J=9.830e-06)
        S['200  x  100']['5.0'] = RHS(1.000e-01, 2.000e-01, 5.000e-03, 7.500e-03, area=2.870e-03, Iyy=1.500e-05, Izz=5.050e-06, J=1.200e-05)
        S['200  x  100']['5.6'] = RHS(1.000e-01, 2.000e-01, 5.600e-03, 8.400e-03, area=3.200e-03, Iyy=1.650e-05, Izz=5.560e-06, J=1.330e-05)
        S['200  x  100']['6.3'] = RHS(1.000e-01, 2.000e-01, 6.300e-03, 9.450e-03, area=3.580e-03, Iyy=1.830e-05, Izz=6.130e-06, J=1.480e-05)
        S['200  x  100']['7.1'] = RHS(1.000e-01, 2.000e-01, 7.100e-03, 1.065e-02, area=4.000e-03, Iyy=2.020e-05, Izz=6.740e-06, J=1.630e-05)
        S['200  x  100']['8.0'] = RHS(1.000e-01, 2.000e-01, 8.000e-03, 1.200e-02, area=4.480e-03, Iyy=2.230e-05, Izz=7.390e-06, J=1.800e-05)
        S['200  x  100']['8.8'] = RHS(1.000e-01, 2.000e-01, 8.800e-03, 1.320e-02, area=4.890e-03, Iyy=2.410e-05, Izz=7.930e-06, J=1.950e-05)
        S['200  x  100']['10.0'] = RHS(1.000e-01, 2.000e-01, 1.000e-02, 1.500e-02, area=5.490e-03, Iyy=2.660e-05, Izz=8.690e-06, J=2.160e-05)
        S['200  x  100']['11.0'] = RHS(1.000e-01, 2.000e-01, 1.100e-02, 1.650e-02, area=5.990e-03, Iyy=2.860e-05, Izz=9.260e-06, J=2.320e-05)
        S['200  x  100']['12.5'] = RHS(1.000e-01, 2.000e-01, 1.250e-02, 1.875e-02, area=6.710e-03, Iyy=3.140e-05, Izz=1.000e-05, J=2.540e-05)
        S['200  x  100']['14.2'] = RHS(1.000e-01, 2.000e-01, 1.420e-02, 2.130e-02, area=7.500e-03, Iyy=3.420e-05, Izz=1.080e-05, J=2.770e-05)
        S['200  x  100']['16.0'] = RHS(1.000e-01, 2.000e-01, 1.600e-02, 2.400e-02, area=8.300e-03, Iyy=3.680e-05, Izz=1.150e-05, J=2.980e-05)
        S['200  x  120'] = {}
        S['200  x  120']['5.0'] = RHS(1.200e-01, 2.000e-01, 5.000e-03, 7.500e-03, area=3.070e-03, Iyy=1.680e-05, Izz=7.620e-06, J=1.650e-05)
        S['200  x  120']['5.6'] = RHS(1.200e-01, 2.000e-01, 5.600e-03, 8.400e-03, area=3.420e-03, Iyy=1.860e-05, Izz=8.410e-06, J=1.830e-05)
        S['200  x  120']['6.3'] = RHS(1.200e-01, 2.000e-01, 6.300e-03, 9.450e-03, area=3.830e-03, Iyy=2.060e-05, Izz=9.290e-06, J=2.030e-05)
        S['200  x  120']['7.1'] = RHS(1.200e-01, 2.000e-01, 7.100e-03, 1.065e-02, area=4.290e-03, Iyy=2.290e-05, Izz=1.020e-05, J=2.250e-05)
        S['200  x  120']['8.0'] = RHS(1.200e-01, 2.000e-01, 8.000e-03, 1.200e-02, area=4.800e-03, Iyy=2.530e-05, Izz=1.130e-05, J=2.500e-05)
        S['200  x  120']['8.8'] = RHS(1.200e-01, 2.000e-01, 8.800e-03, 1.320e-02, area=5.240e-03, Iyy=2.730e-05, Izz=1.220e-05, J=2.700e-05)
        S['200  x  120']['10.0'] = RHS(1.200e-01, 2.000e-01, 1.000e-02, 1.500e-02, area=5.890e-03, Iyy=3.030e-05, Izz=1.340e-05, J=3.000e-05)
        S['200  x  120']['11.0'] = RHS(1.200e-01, 2.000e-01, 1.100e-02, 1.650e-02, area=6.430e-03, Iyy=3.260e-05, Izz=1.430e-05, J=3.240e-05)
        S['200  x  120']['12.5'] = RHS(1.200e-01, 2.000e-01, 1.250e-02, 1.875e-02, area=7.210e-03, Iyy=3.580e-05, Izz=1.560e-05, J=3.570e-05)
        S['200  x  120']['14.2'] = RHS(1.200e-01, 2.000e-01, 1.420e-02, 2.130e-02, area=8.070e-03, Iyy=3.910e-05, Izz=1.690e-05, J=3.920e-05)
        S['200  x  120']['16.0'] = RHS(1.200e-01, 2.000e-01, 1.600e-02, 2.400e-02, area=8.940e-03, Iyy=4.220e-05, Izz=1.810e-05, J=4.250e-05)
        S['200  x  150'] = {}
        S['200  x  150']['5.0'] = RHS(1.500e-01, 2.000e-01, 5.000e-03, 7.500e-03, area=3.370e-03, Iyy=1.970e-05, Izz=1.260e-05, J=2.390e-05)
        S['200  x  150']['5.6'] = RHS(1.500e-01, 2.000e-01, 5.600e-03, 8.400e-03, area=3.760e-03, Iyy=2.180e-05, Izz=1.400e-05, J=2.650e-05)
        S['200  x  150']['6.3'] = RHS(1.500e-01, 2.000e-01, 6.300e-03, 9.450e-03, area=4.210e-03, Iyy=2.420e-05, Izz=1.550e-05, J=2.950e-05)
        S['200  x  150']['7.1'] = RHS(1.500e-01, 2.000e-01, 7.100e-03, 1.065e-02, area=4.710e-03, Iyy=2.680e-05, Izz=1.720e-05, J=3.280e-05)
        S['200  x  150']['8.0'] = RHS(1.500e-01, 2.000e-01, 8.000e-03, 1.200e-02, area=5.280e-03, Iyy=2.970e-05, Izz=1.890e-05, J=3.640e-05)
        S['200  x  150']['8.8'] = RHS(1.500e-01, 2.000e-01, 8.800e-03, 1.320e-02, area=5.770e-03, Iyy=3.220e-05, Izz=2.050e-05, J=3.960e-05)
        S['200  x  150']['10.0'] = RHS(1.500e-01, 2.000e-01, 1.000e-02, 1.500e-02, area=6.490e-03, Iyy=3.570e-05, Izz=2.260e-05, J=4.410e-05)
        S['200  x  150']['11.0'] = RHS(1.500e-01, 2.000e-01, 1.100e-02, 1.650e-02, area=7.090e-03, Iyy=3.840e-05, Izz=2.440e-05, J=4.770e-05)
        S['200  x  150']['12.5'] = RHS(1.500e-01, 2.000e-01, 1.250e-02, 1.875e-02, area=7.960e-03, Iyy=4.240e-05, Izz=2.670e-05, J=5.290e-05)
        S['200  x  150']['14.2'] = RHS(1.500e-01, 2.000e-01, 1.420e-02, 2.130e-02, area=8.920e-03, Iyy=4.640e-05, Izz=2.920e-05, J=5.830e-05)
        S['200  x  150']['16.0'] = RHS(1.500e-01, 2.000e-01, 1.600e-02, 2.400e-02, area=9.900e-03, Iyy=5.040e-05, Izz=3.150e-05, J=6.370e-05)
        S['220  x  120'] = {}
        S['220  x  120']['7.1'] = RHS(1.200e-01, 2.200e-01, 7.100e-03, 1.065e-02, area=4.570e-03, Iyy=2.900e-05, Izz=1.120e-05, J=2.570e-05)
        S['220  x  120']['8.0'] = RHS(1.200e-01, 2.200e-01, 8.000e-03, 1.200e-02, area=5.120e-03, Iyy=3.200e-05, Izz=1.230e-05, J=2.850e-05)
        S['220  x  120']['8.8'] = RHS(1.200e-01, 2.200e-01, 8.800e-03, 1.320e-02, area=5.590e-03, Iyy=3.470e-05, Izz=1.320e-05, J=3.090e-05)
        S['220  x  120']['10.0'] = RHS(1.200e-01, 2.200e-01, 1.000e-02, 1.500e-02, area=6.290e-03, Iyy=3.840e-05, Izz=1.460e-05, J=3.430e-05)
        S['220  x  120']['11.0'] = RHS(1.200e-01, 2.200e-01, 1.100e-02, 1.650e-02, area=6.870e-03, Iyy=4.140e-05, Izz=1.560e-05, J=3.700e-05)
        S['220  x  120']['12.5'] = RHS(1.200e-01, 2.200e-01, 1.250e-02, 1.875e-02, area=7.710e-03, Iyy=4.560e-05, Izz=1.710e-05, J=4.090e-05)
        S['220  x  120']['14.2'] = RHS(1.200e-01, 2.200e-01, 1.420e-02, 2.130e-02, area=8.630e-03, Iyy=5.000e-05, Izz=1.850e-05, J=4.490e-05)
        S['220  x  120']['16.0'] = RHS(1.200e-01, 2.200e-01, 1.600e-02, 2.400e-02, area=9.580e-03, Iyy=5.410e-05, Izz=1.990e-05, J=4.870e-05)
        S['250  x  100'] = {}
        S['250  x  100']['5.0'] = RHS(1.000e-01, 2.500e-01, 5.000e-03, 7.500e-03, area=3.370e-03, Iyy=2.610e-05, Izz=6.180e-06, J=1.620e-05)
        S['250  x  100']['5.6'] = RHS(1.000e-01, 2.500e-01, 5.600e-03, 8.400e-03, area=3.760e-03, Iyy=2.890e-05, Izz=6.810e-06, J=1.790e-05)
        S['250  x  100']['6.3'] = RHS(1.000e-01, 2.500e-01, 6.300e-03, 9.450e-03, area=4.210e-03, Iyy=3.210e-05, Izz=7.510e-06, J=1.980e-05)
        S['250  x  100']['7.1'] = RHS(1.000e-01, 2.500e-01, 7.100e-03, 1.065e-02, area=4.710e-03, Iyy=3.560e-05, Izz=8.270e-06, J=2.200e-05)
        S['250  x  100']['8.0'] = RHS(1.000e-01, 2.500e-01, 8.000e-03, 1.200e-02, area=5.280e-03, Iyy=3.940e-05, Izz=9.090e-06, J=2.430e-05)
        S['250  x  100']['8.8'] = RHS(1.000e-01, 2.500e-01, 8.800e-03, 1.320e-02, area=5.770e-03, Iyy=4.270e-05, Izz=9.770e-06, J=2.630e-05)
        S['250  x  100']['10.0'] = RHS(1.000e-01, 2.500e-01, 1.000e-02, 1.500e-02, area=6.490e-03, Iyy=4.730e-05, Izz=1.070e-05, J=2.910e-05)
        S['250  x  100']['11.0'] = RHS(1.000e-01, 2.500e-01, 1.100e-02, 1.650e-02, area=7.090e-03, Iyy=5.100e-05, Izz=1.140e-05, J=3.130e-05)
        S['250  x  100']['12.5'] = RHS(1.000e-01, 2.500e-01, 1.250e-02, 1.875e-02, area=7.960e-03, Iyy=5.620e-05, Izz=1.240e-05, J=3.440e-05)
        S['250  x  100']['14.2'] = RHS(1.000e-01, 2.500e-01, 1.420e-02, 2.130e-02, area=8.920e-03, Iyy=6.160e-05, Izz=1.340e-05, J=3.750e-05)
        S['250  x  100']['16.0'] = RHS(1.000e-01, 2.500e-01, 1.600e-02, 2.400e-02, area=9.900e-03, Iyy=6.690e-05, Izz=1.430e-05, J=4.050e-05)
        S['250  x  150'] = {}
        S['250  x  150']['5.0'] = RHS(1.500e-01, 2.500e-01, 5.000e-03, 7.500e-03, area=3.870e-03, Iyy=3.360e-05, Izz=1.530e-05, J=3.280e-05)
        S['250  x  150']['5.6'] = RHS(1.500e-01, 2.500e-01, 5.600e-03, 8.400e-03, area=4.320e-03, Iyy=3.730e-05, Izz=1.690e-05, J=3.640e-05)
        S['250  x  150']['6.3'] = RHS(1.500e-01, 2.500e-01, 6.300e-03, 9.450e-03, area=4.840e-03, Iyy=4.140e-05, Izz=1.870e-05, J=4.050e-05)
        S['250  x  150']['7.1'] = RHS(1.500e-01, 2.500e-01, 7.100e-03, 1.065e-02, area=5.420e-03, Iyy=4.610e-05, Izz=2.080e-05, J=4.520e-05)
        S['250  x  150']['8.0'] = RHS(1.500e-01, 2.500e-01, 8.000e-03, 1.200e-02, area=6.080e-03, Iyy=5.110e-05, Izz=2.300e-05, J=5.020e-05)
        S['250  x  150']['8.8'] = RHS(1.500e-01, 2.500e-01, 8.800e-03, 1.320e-02, area=6.650e-03, Iyy=5.550e-05, Izz=2.490e-05, J=5.460e-05)
        S['250  x  150']['10.0'] = RHS(1.500e-01, 2.500e-01, 1.000e-02, 1.500e-02, area=7.490e-03, Iyy=6.170e-05, Izz=2.760e-05, J=6.090e-05)
        S['250  x  150']['11.0'] = RHS(1.500e-01, 2.500e-01, 1.100e-02, 1.650e-02, area=8.190e-03, Iyy=6.670e-05, Izz=2.970e-05, J=6.600e-05)
        S['250  x  150']['12.5'] = RHS(1.500e-01, 2.500e-01, 1.250e-02, 1.875e-02, area=9.210e-03, Iyy=7.390e-05, Izz=3.260e-05, J=7.330e-05)
        S['250  x  150']['14.2'] = RHS(1.500e-01, 2.500e-01, 1.420e-02, 2.130e-02, area=1.030e-02, Iyy=8.140e-05, Izz=3.580e-05, J=8.100e-05)
        S['250  x  150']['16.0'] = RHS(1.500e-01, 2.500e-01, 1.600e-02, 2.400e-02, area=1.150e-02, Iyy=8.880e-05, Izz=3.870e-05, J=8.870e-05)
        S['260  x  140'] = {}
        S['260  x  140']['5.0'] = RHS(1.400e-01, 2.600e-01, 5.000e-03, 7.500e-03, area=3.870e-03, Iyy=3.530e-05, Izz=1.350e-05, J=3.080e-05)
        S['260  x  140']['5.6'] = RHS(1.400e-01, 2.600e-01, 5.600e-03, 8.400e-03, area=4.320e-03, Iyy=3.920e-05, Izz=1.500e-05, J=3.420e-05)
        S['260  x  140']['6.3'] = RHS(1.400e-01, 2.600e-01, 6.300e-03, 9.450e-03, area=4.840e-03, Iyy=4.360e-05, Izz=1.660e-05, J=3.800e-05)
        S['260  x  140']['7.1'] = RHS(1.400e-01, 2.600e-01, 7.100e-03, 1.065e-02, area=5.420e-03, Iyy=4.840e-05, Izz=1.840e-05, J=4.230e-05)
        S['260  x  140']['8.0'] = RHS(1.400e-01, 2.600e-01, 8.000e-03, 1.200e-02, area=6.080e-03, Iyy=5.370e-05, Izz=2.030e-05, J=4.700e-05)
        S['260  x  140']['8.8'] = RHS(1.400e-01, 2.600e-01, 8.800e-03, 1.320e-02, area=6.650e-03, Iyy=5.830e-05, Izz=2.200e-05, J=5.110e-05)
        S['260  x  140']['10.0'] = RHS(1.400e-01, 2.600e-01, 1.000e-02, 1.500e-02, area=7.490e-03, Iyy=6.490e-05, Izz=2.430e-05, J=5.700e-05)
        S['260  x  140']['11.0'] = RHS(1.400e-01, 2.600e-01, 1.100e-02, 1.650e-02, area=8.190e-03, Iyy=7.020e-05, Izz=2.620e-05, J=6.170e-05)
        S['260  x  140']['12.5'] = RHS(1.400e-01, 2.600e-01, 1.250e-02, 1.875e-02, area=9.210e-03, Iyy=7.770e-05, Izz=2.880e-05, J=6.840e-05)
        S['260  x  140']['14.2'] = RHS(1.400e-01, 2.600e-01, 1.420e-02, 2.130e-02, area=1.030e-02, Iyy=8.560e-05, Izz=3.140e-05, J=7.560e-05)
        S['260  x  140']['16.0'] = RHS(1.400e-01, 2.600e-01, 1.600e-02, 2.400e-02, area=1.150e-02, Iyy=9.340e-05, Izz=3.400e-05, J=8.260e-05)
        S['260  x  180'] = {}
        S['260  x  180']['8.0'] = RHS(1.800e-01, 2.600e-01, 8.000e-03, 1.200e-02, area=6.720e-03, Iyy=6.390e-05, Izz=3.610e-05, J=7.220e-05)
        S['260  x  180']['8.8'] = RHS(1.800e-01, 2.600e-01, 8.800e-03, 1.320e-02, area=7.350e-03, Iyy=6.940e-05, Izz=3.910e-05, J=7.860e-05)
        S['260  x  180']['10.0'] = RHS(1.800e-01, 2.600e-01, 1.000e-02, 1.500e-02, area=8.290e-03, Iyy=7.740e-05, Izz=4.350e-05, J=8.800e-05)
        S['260  x  180']['11.0'] = RHS(1.800e-01, 2.600e-01, 1.100e-02, 1.650e-02, area=9.070e-03, Iyy=8.380e-05, Izz=4.700e-05, J=9.550e-05)
        S['260  x  180']['12.5'] = RHS(1.800e-01, 2.600e-01, 1.250e-02, 1.875e-02, area=1.020e-02, Iyy=9.300e-05, Izz=5.200e-05, J=1.060e-04)
        S['260  x  180']['14.2'] = RHS(1.800e-01, 2.600e-01, 1.420e-02, 2.130e-02, area=1.150e-02, Iyy=1.030e-04, Izz=5.720e-05, J=1.180e-04)
        S['260  x  180']['16.0'] = RHS(1.800e-01, 2.600e-01, 1.600e-02, 2.400e-02, area=1.280e-02, Iyy=1.120e-04, Izz=6.230e-05, J=1.300e-04)
        S['300  x  100'] = {}
        S['300  x  100']['5.0'] = RHS(1.000e-01, 3.000e-01, 5.000e-03, 7.500e-03, area=3.870e-03, Iyy=4.150e-05, Izz=7.310e-06, J=2.040e-05)
        S['300  x  100']['5.6'] = RHS(1.000e-01, 3.000e-01, 5.600e-03, 8.400e-03, area=4.320e-03, Iyy=4.600e-05, Izz=8.060e-06, J=2.260e-05)
        S['300  x  100']['6.3'] = RHS(1.000e-01, 3.000e-01, 6.300e-03, 9.450e-03, area=4.840e-03, Iyy=5.110e-05, Izz=8.900e-06, J=2.500e-05)
        S['300  x  100']['7.1'] = RHS(1.000e-01, 3.000e-01, 7.100e-03, 1.065e-02, area=5.420e-03, Iyy=5.680e-05, Izz=9.810e-06, J=2.780e-05)
        S['300  x  100']['8.0'] = RHS(1.000e-01, 3.000e-01, 8.000e-03, 1.200e-02, area=6.080e-03, Iyy=6.300e-05, Izz=1.080e-05, J=3.070e-05)
        S['300  x  100']['8.8'] = RHS(1.000e-01, 3.000e-01, 8.800e-03, 1.320e-02, area=6.650e-03, Iyy=6.840e-05, Izz=1.160e-05, J=3.320e-05)
        S['300  x  100']['10.0'] = RHS(1.000e-01, 3.000e-01, 1.000e-02, 1.500e-02, area=7.490e-03, Iyy=7.610e-05, Izz=1.280e-05, J=3.680e-05)
        S['300  x  100']['11.0'] = RHS(1.000e-01, 3.000e-01, 1.100e-02, 1.650e-02, area=8.190e-03, Iyy=8.230e-05, Izz=1.360e-05, J=3.960e-05)
        S['300  x  100']['12.5'] = RHS(1.000e-01, 3.000e-01, 1.250e-02, 1.875e-02, area=9.210e-03, Iyy=9.100e-05, Izz=1.490e-05, J=4.350e-05)
        S['300  x  100']['14.2'] = RHS(1.000e-01, 3.000e-01, 1.420e-02, 2.130e-02, area=1.030e-02, Iyy=1.000e-04, Izz=1.610e-05, J=4.760e-05)
        S['300  x  100']['16.0'] = RHS(1.000e-01, 3.000e-01, 1.600e-02, 2.400e-02, area=1.150e-02, Iyy=1.090e-04, Izz=1.720e-05, J=5.140e-05)
        S['300  x  200'] = {}
        S['300  x  200']['5.0'] = RHS(2.000e-01, 3.000e-01, 5.000e-03, 7.500e-03, area=4.870e-03, Iyy=6.320e-05, Izz=3.400e-05, J=6.820e-05)
        S['300  x  200']['5.6'] = RHS(2.000e-01, 3.000e-01, 5.600e-03, 8.400e-03, area=5.440e-03, Iyy=7.020e-05, Izz=3.770e-05, J=7.590e-05)
        S['300  x  200']['6.3'] = RHS(2.000e-01, 3.000e-01, 6.300e-03, 9.450e-03, area=6.100e-03, Iyy=7.830e-05, Izz=4.190e-05, J=8.480e-05)
        S['300  x  200']['7.1'] = RHS(2.000e-01, 3.000e-01, 7.100e-03, 1.065e-02, area=6.840e-03, Iyy=8.730e-05, Izz=4.670e-05, J=9.470e-05)
        S['300  x  200']['8.0'] = RHS(2.000e-01, 3.000e-01, 8.000e-03, 1.200e-02, area=7.680e-03, Iyy=9.720e-05, Izz=5.180e-05, J=1.060e-04)
        S['300  x  200']['8.8'] = RHS(2.000e-01, 3.000e-01, 8.800e-03, 1.320e-02, area=8.410e-03, Iyy=1.060e-04, Izz=5.630e-05, J=1.150e-04)
        S['300  x  200']['10.0'] = RHS(2.000e-01, 3.000e-01, 1.000e-02, 1.500e-02, area=9.490e-03, Iyy=1.180e-04, Izz=6.280e-05, J=1.290e-04)
        S['300  x  200']['11.0'] = RHS(2.000e-01, 3.000e-01, 1.100e-02, 1.650e-02, area=1.040e-02, Iyy=1.280e-04, Izz=6.800e-05, J=1.400e-04)
        S['300  x  200']['12.5'] = RHS(2.000e-01, 3.000e-01, 1.250e-02, 1.875e-02, area=1.170e-02, Iyy=1.430e-04, Izz=7.540e-05, J=1.570e-04)
        S['300  x  200']['14.2'] = RHS(2.000e-01, 3.000e-01, 1.420e-02, 2.130e-02, area=1.320e-02, Iyy=1.580e-04, Izz=8.330e-05, J=1.750e-04)
        S['300  x  200']['16.0'] = RHS(2.000e-01, 3.000e-01, 1.600e-02, 2.400e-02, area=1.470e-02, Iyy=1.740e-04, Izz=9.110e-05, J=1.930e-04)
        S['300  x  250'] = {}
        S['300  x  250']['6.3'] = RHS(2.500e-01, 3.000e-01, 6.300e-03, 9.450e-03, area=6.730e-03, Iyy=9.190e-05, Izz=6.950e-05, J=1.220e-04)
        S['300  x  250']['7.1'] = RHS(2.500e-01, 3.000e-01, 7.100e-03, 1.065e-02, area=7.550e-03, Iyy=1.030e-04, Izz=7.750e-05, J=1.360e-04)
        S['300  x  250']['8.0'] = RHS(2.500e-01, 3.000e-01, 8.000e-03, 1.200e-02, area=8.480e-03, Iyy=1.140e-04, Izz=8.630e-05, J=1.520e-04)
        S['300  x  250']['8.8'] = RHS(2.500e-01, 3.000e-01, 8.800e-03, 1.320e-02, area=9.290e-03, Iyy=1.240e-04, Izz=9.390e-05, J=1.660e-04)
        S['300  x  250']['10.0'] = RHS(2.500e-01, 3.000e-01, 1.000e-02, 1.500e-02, area=1.050e-02, Iyy=1.390e-04, Izz=1.050e-04, J=1.860e-04)
        S['300  x  250']['11.0'] = RHS(2.500e-01, 3.000e-01, 1.100e-02, 1.650e-02, area=1.150e-02, Iyy=1.510e-04, Izz=1.140e-04, J=2.030e-04)
        S['300  x  250']['12.5'] = RHS(2.500e-01, 3.000e-01, 1.250e-02, 1.875e-02, area=1.300e-02, Iyy=1.690e-04, Izz=1.270e-04, J=2.270e-04)
        S['300  x  250']['14.2'] = RHS(2.500e-01, 3.000e-01, 1.420e-02, 2.130e-02, area=1.460e-02, Iyy=1.870e-04, Izz=1.410e-04, J=2.540e-04)
        S['300  x  250']['16.0'] = RHS(2.500e-01, 3.000e-01, 1.600e-02, 2.400e-02, area=1.630e-02, Iyy=2.060e-04, Izz=1.550e-04, J=2.810e-04)
        S['350  x  150'] = {}
        S['350  x  150']['5.0'] = RHS(1.500e-01, 3.500e-01, 5.000e-03, 7.500e-03, area=4.870e-03, Iyy=7.660e-05, Izz=2.050e-05, J=5.160e-05)
        S['350  x  150']['5.6'] = RHS(1.500e-01, 3.500e-01, 5.600e-03, 8.400e-03, area=5.440e-03, Iyy=8.510e-05, Izz=2.270e-05, J=5.730e-05)
        S['350  x  150']['6.3'] = RHS(1.500e-01, 3.500e-01, 6.300e-03, 9.450e-03, area=6.100e-03, Iyy=9.480e-05, Izz=2.520e-05, J=6.390e-05)
        S['350  x  150']['7.1'] = RHS(1.500e-01, 3.500e-01, 7.100e-03, 1.065e-02, area=6.840e-03, Iyy=1.060e-04, Izz=2.800e-05, J=7.120e-05)
        S['350  x  150']['8.0'] = RHS(1.500e-01, 3.500e-01, 8.000e-03, 1.200e-02, area=7.680e-03, Iyy=1.180e-04, Izz=3.100e-05, J=7.930e-05)
        S['350  x  150']['8.8'] = RHS(1.500e-01, 3.500e-01, 8.800e-03, 1.320e-02, area=8.410e-03, Iyy=1.280e-04, Izz=3.360e-05, J=8.620e-05)
        S['350  x  150']['10.0'] = RHS(1.500e-01, 3.500e-01, 1.000e-02, 1.500e-02, area=9.490e-03, Iyy=1.430e-04, Izz=3.740e-05, J=9.630e-05)
        S['350  x  150']['11.0'] = RHS(1.500e-01, 3.500e-01, 1.100e-02, 1.650e-02, area=1.040e-02, Iyy=1.550e-04, Izz=4.030e-05, J=1.040e-04)
        S['350  x  150']['12.5'] = RHS(1.500e-01, 3.500e-01, 1.250e-02, 1.875e-02, area=1.170e-02, Iyy=1.730e-04, Izz=4.450e-05, J=1.160e-04)
        S['350  x  150']['14.2'] = RHS(1.500e-01, 3.500e-01, 1.420e-02, 2.130e-02, area=1.320e-02, Iyy=1.920e-04, Izz=4.890e-05, J=1.290e-04)
        S['350  x  150']['16.0'] = RHS(1.500e-01, 3.500e-01, 1.600e-02, 2.400e-02, area=1.470e-02, Iyy=2.110e-04, Izz=5.320e-05, J=1.410e-04)
        S['350  x  200'] = {}
        S['350  x  200']['6.3'] = RHS(2.000e-01, 3.500e-01, 6.300e-03, 9.450e-03, area=6.730e-03, Iyy=1.130e-04, Izz=4.780e-05, J=1.050e-04)
        S['350  x  200']['7.1'] = RHS(2.000e-01, 3.500e-01, 7.100e-03, 1.065e-02, area=7.550e-03, Iyy=1.270e-04, Izz=5.330e-05, J=1.180e-04)
        S['350  x  200']['8.0'] = RHS(2.000e-01, 3.500e-01, 8.000e-03, 1.200e-02, area=8.480e-03, Iyy=1.410e-04, Izz=5.920e-05, J=1.310e-04)
        S['350  x  200']['8.8'] = RHS(2.000e-01, 3.500e-01, 8.800e-03, 1.320e-02, area=9.290e-03, Iyy=1.540e-04, Izz=6.440e-05, J=1.430e-04)
        S['350  x  200']['10.0'] = RHS(2.000e-01, 3.500e-01, 1.000e-02, 1.500e-02, area=1.050e-02, Iyy=1.720e-04, Izz=7.180e-05, J=1.600e-04)
        S['350  x  200']['11.0'] = RHS(2.000e-01, 3.500e-01, 1.100e-02, 1.650e-02, area=1.150e-02, Iyy=1.870e-04, Izz=7.780e-05, J=1.750e-04)
        S['350  x  200']['12.5'] = RHS(2.000e-01, 3.500e-01, 1.250e-02, 1.875e-02, area=1.300e-02, Iyy=2.090e-04, Izz=8.640e-05, J=1.950e-04)
        S['350  x  200']['14.2'] = RHS(2.000e-01, 3.500e-01, 1.420e-02, 2.130e-02, area=1.460e-02, Iyy=2.320e-04, Izz=9.560e-05, J=2.170e-04)
        S['350  x  200']['16.0'] = RHS(2.000e-01, 3.500e-01, 1.600e-02, 2.400e-02, area=1.630e-02, Iyy=2.550e-04, Izz=1.050e-04, J=2.400e-04)
        S['350  x  250'] = {}
        S['350  x  250']['6.3'] = RHS(2.500e-01, 3.500e-01, 6.300e-03, 9.450e-03, area=7.360e-03, Iyy=1.320e-04, Izz=7.880e-05, J=1.520e-04)
        S['350  x  250']['7.1'] = RHS(2.500e-01, 3.500e-01, 7.100e-03, 1.065e-02, area=8.260e-03, Iyy=1.470e-04, Izz=8.800e-05, J=1.700e-04)
        S['350  x  250']['8.0'] = RHS(2.500e-01, 3.500e-01, 8.000e-03, 1.200e-02, area=9.280e-03, Iyy=1.640e-04, Izz=9.800e-05, J=1.900e-04)
        S['350  x  250']['8.8'] = RHS(2.500e-01, 3.500e-01, 8.800e-03, 1.320e-02, area=1.020e-02, Iyy=1.790e-04, Izz=1.070e-04, J=2.080e-04)
        S['350  x  250']['10.0'] = RHS(2.500e-01, 3.500e-01, 1.000e-02, 1.500e-02, area=1.150e-02, Iyy=2.010e-04, Izz=1.190e-04, J=2.340e-04)
        S['350  x  250']['11.0'] = RHS(2.500e-01, 3.500e-01, 1.100e-02, 1.650e-02, area=1.260e-02, Iyy=2.190e-04, Izz=1.300e-04, J=2.550e-04)
        S['350  x  250']['12.5'] = RHS(2.500e-01, 3.500e-01, 1.250e-02, 1.875e-02, area=1.420e-02, Iyy=2.440e-04, Izz=1.440e-04, J=2.850e-04)
        S['350  x  250']['14.2'] = RHS(2.500e-01, 3.500e-01, 1.420e-02, 2.130e-02, area=1.600e-02, Iyy=2.720e-04, Izz=1.600e-04, J=3.190e-04)
        S['350  x  250']['16.0'] = RHS(2.500e-01, 3.500e-01, 1.600e-02, 2.400e-02, area=1.790e-02, Iyy=3.000e-04, Izz=1.770e-04, J=3.530e-04)
        S['400  x  150'] = {}
        S['400  x  150']['6.3'] = RHS(1.500e-01, 4.000e-01, 6.300e-03, 9.450e-03, area=6.730e-03, Iyy=1.330e-04, Izz=2.850e-05, J=7.600e-05)
        S['400  x  150']['7.1'] = RHS(1.500e-01, 4.000e-01, 7.100e-03, 1.065e-02, area=7.550e-03, Iyy=1.480e-04, Izz=3.170e-05, J=8.470e-05)
        S['400  x  150']['8.0'] = RHS(1.500e-01, 4.000e-01, 8.000e-03, 1.200e-02, area=8.480e-03, Iyy=1.650e-04, Izz=3.510e-05, J=9.420e-05)
        S['400  x  150']['8.8'] = RHS(1.500e-01, 4.000e-01, 8.800e-03, 1.320e-02, area=9.290e-03, Iyy=1.800e-04, Izz=3.800e-05, J=1.030e-04)
        S['400  x  150']['10.0'] = RHS(1.500e-01, 4.000e-01, 1.000e-02, 1.500e-02, area=1.050e-02, Iyy=2.010e-04, Izz=4.230e-05, J=1.150e-04)
        S['400  x  150']['11.0'] = RHS(1.500e-01, 4.000e-01, 1.100e-02, 1.650e-02, area=1.150e-02, Iyy=2.180e-04, Izz=4.560e-05, J=1.240e-04)
        S['400  x  150']['12.5'] = RHS(1.500e-01, 4.000e-01, 1.250e-02, 1.875e-02, area=1.300e-02, Iyy=2.440e-04, Izz=5.040e-05, J=1.380e-04)
        S['400  x  150']['14.2'] = RHS(1.500e-01, 4.000e-01, 1.420e-02, 2.130e-02, area=1.460e-02, Iyy=2.710e-04, Izz=5.550e-05, J=1.530e-04)
        S['400  x  150']['16.0'] = RHS(1.500e-01, 4.000e-01, 1.600e-02, 2.400e-02, area=1.630e-02, Iyy=2.980e-04, Izz=6.040e-05, J=1.680e-04)
        S['400  x  200'] = {}
        S['400  x  200']['6.3'] = RHS(2.000e-01, 4.000e-01, 6.300e-03, 9.450e-03, area=7.360e-03, Iyy=1.570e-04, Izz=5.380e-05, J=1.260e-04)
        S['400  x  200']['7.1'] = RHS(2.000e-01, 4.000e-01, 7.100e-03, 1.065e-02, area=8.260e-03, Iyy=1.750e-04, Izz=5.990e-05, J=1.410e-04)
        S['400  x  200']['8.0'] = RHS(2.000e-01, 4.000e-01, 8.000e-03, 1.200e-02, area=9.280e-03, Iyy=1.960e-04, Izz=6.660e-05, J=1.570e-04)
        S['400  x  200']['8.8'] = RHS(2.000e-01, 4.000e-01, 8.800e-03, 1.320e-02, area=1.020e-02, Iyy=2.130e-04, Izz=7.240e-05, J=1.720e-04)
        S['400  x  200']['10.0'] = RHS(2.000e-01, 4.000e-01, 1.000e-02, 1.500e-02, area=1.150e-02, Iyy=2.390e-04, Izz=8.080e-05, J=1.930e-04)
        S['400  x  200']['11.0'] = RHS(2.000e-01, 4.000e-01, 1.100e-02, 1.650e-02, area=1.260e-02, Iyy=2.600e-04, Izz=8.760e-05, J=2.100e-04)
        S['400  x  200']['12.5'] = RHS(2.000e-01, 4.000e-01, 1.250e-02, 1.875e-02, area=1.420e-02, Iyy=2.910e-04, Izz=9.740e-05, J=2.340e-04)
        S['400  x  200']['14.2'] = RHS(2.000e-01, 4.000e-01, 1.420e-02, 2.130e-02, area=1.600e-02, Iyy=3.240e-04, Izz=1.080e-04, J=2.610e-04)
        S['400  x  200']['16.0'] = RHS(2.000e-01, 4.000e-01, 1.600e-02, 2.400e-02, area=1.790e-02, Iyy=3.570e-04, Izz=1.180e-04, J=2.890e-04)
        S['400  x  300'] = {}
        S['400  x  300']['8.0'] = RHS(3.000e-01, 4.000e-01, 8.000e-03, 1.200e-02, area=1.090e-02, Iyy=2.570e-04, Izz=1.650e-04, J=3.100e-04)
        S['400  x  300']['8.8'] = RHS(3.000e-01, 4.000e-01, 8.800e-03, 1.320e-02, area=1.190e-02, Iyy=2.810e-04, Izz=1.800e-04, J=3.390e-04)
        S['400  x  300']['10.0'] = RHS(3.000e-01, 4.000e-01, 1.000e-02, 1.500e-02, area=1.350e-02, Iyy=3.150e-04, Izz=2.020e-04, J=3.820e-04)
        S['400  x  300']['11.0'] = RHS(3.000e-01, 4.000e-01, 1.100e-02, 1.650e-02, area=1.480e-02, Iyy=3.430e-04, Izz=2.200e-04, J=4.170e-04)
        S['400  x  300']['12.5'] = RHS(3.000e-01, 4.000e-01, 1.250e-02, 1.875e-02, area=1.670e-02, Iyy=3.850e-04, Izz=2.460e-04, J=4.680e-04)
        S['400  x  300']['14.2'] = RHS(3.000e-01, 4.000e-01, 1.420e-02, 2.130e-02, area=1.890e-02, Iyy=4.300e-04, Izz=2.740e-04, J=5.250e-04)
        S['400  x  300']['16.0'] = RHS(3.000e-01, 4.000e-01, 1.600e-02, 2.400e-02, area=2.110e-02, Iyy=4.750e-04, Izz=3.030e-04, J=5.830e-04)
        S['450  x  150'] = {}
        S['450  x  150']['8.0'] = RHS(1.500e-01, 4.500e-01, 8.000e-03, 1.200e-02, area=9.280e-03, Iyy=2.230e-04, Izz=3.910e-05, J=1.090e-04)
        S['450  x  150']['8.8'] = RHS(1.500e-01, 4.500e-01, 8.800e-03, 1.320e-02, area=1.020e-02, Iyy=2.430e-04, Izz=4.240e-05, J=1.190e-04)
        S['450  x  150']['10.0'] = RHS(1.500e-01, 4.500e-01, 1.000e-02, 1.500e-02, area=1.150e-02, Iyy=2.720e-04, Izz=4.720e-05, J=1.330e-04)
        S['450  x  150']['11.0'] = RHS(1.500e-01, 4.500e-01, 1.100e-02, 1.650e-02, area=1.260e-02, Iyy=2.960e-04, Izz=5.100e-05, J=1.440e-04)
        S['450  x  150']['12.5'] = RHS(1.500e-01, 4.500e-01, 1.250e-02, 1.875e-02, area=1.420e-02, Iyy=3.310e-04, Izz=5.640e-05, J=1.610e-04)
        S['450  x  150']['14.2'] = RHS(1.500e-01, 4.500e-01, 1.420e-02, 2.130e-02, area=1.600e-02, Iyy=3.680e-04, Izz=6.200e-05, J=1.780e-04)
        S['450  x  150']['16.0'] = RHS(1.500e-01, 4.500e-01, 1.600e-02, 2.400e-02, area=1.790e-02, Iyy=4.060e-04, Izz=6.760e-05, J=1.960e-04)
        S['450  x  250'] = {}
        S['450  x  250']['8.0'] = RHS(2.500e-01, 4.500e-01, 8.000e-03, 1.200e-02, area=1.090e-02, Iyy=3.010e-04, Izz=1.210e-04, J=2.710e-04)
        S['450  x  250']['8.8'] = RHS(2.500e-01, 4.500e-01, 8.800e-03, 1.320e-02, area=1.190e-02, Iyy=3.280e-04, Izz=1.320e-04, J=2.960e-04)
        S['450  x  250']['10.0'] = RHS(2.500e-01, 4.500e-01, 1.000e-02, 1.500e-02, area=1.350e-02, Iyy=3.690e-04, Izz=1.480e-04, J=3.330e-04)
        S['450  x  250']['11.0'] = RHS(2.500e-01, 4.500e-01, 1.100e-02, 1.650e-02, area=1.480e-02, Iyy=4.020e-04, Izz=1.610e-04, J=3.630e-04)
        S['450  x  250']['12.5'] = RHS(2.500e-01, 4.500e-01, 1.250e-02, 1.875e-02, area=1.670e-02, Iyy=4.500e-04, Izz=1.800e-04, J=4.070e-04)
        S['450  x  250']['14.2'] = RHS(2.500e-01, 4.500e-01, 1.420e-02, 2.130e-02, area=1.890e-02, Iyy=5.030e-04, Izz=2.000e-04, J=4.560e-04)
        S['450  x  250']['16.0'] = RHS(2.500e-01, 4.500e-01, 1.600e-02, 2.400e-02, area=2.110e-02, Iyy=5.570e-04, Izz=2.200e-04, J=5.050e-04)
        S['500  x  200'] = {}
        S['500  x  200']['8.0'] = RHS(2.000e-01, 5.000e-01, 8.000e-03, 1.200e-02, area=1.090e-02, Iyy=3.400e-04, Izz=8.140e-05, J=2.110e-04)
        S['500  x  200']['8.8'] = RHS(2.000e-01, 5.000e-01, 8.800e-03, 1.320e-02, area=1.190e-02, Iyy=3.720e-04, Izz=8.850e-05, J=2.300e-04)
        S['500  x  200']['10.0'] = RHS(2.000e-01, 5.000e-01, 1.000e-02, 1.500e-02, area=1.350e-02, Iyy=4.180e-04, Izz=9.890e-05, J=2.590e-04)
        S['500  x  200']['11.0'] = RHS(2.000e-01, 5.000e-01, 1.100e-02, 1.650e-02, area=1.480e-02, Iyy=4.550e-04, Izz=1.070e-04, J=2.820e-04)
        S['500  x  200']['12.5'] = RHS(2.000e-01, 5.000e-01, 1.250e-02, 1.875e-02, area=1.670e-02, Iyy=5.100e-04, Izz=1.190e-04, J=3.150e-04)
        S['500  x  200']['14.2'] = RHS(2.000e-01, 5.000e-01, 1.420e-02, 2.130e-02, area=1.890e-02, Iyy=5.690e-04, Izz=1.320e-04, J=3.520e-04)
        S['500  x  200']['16.0'] = RHS(2.000e-01, 5.000e-01, 1.600e-02, 2.400e-02, area=2.110e-02, Iyy=6.300e-04, Izz=1.450e-04, J=3.890e-04)
        S['500  x  300'] = {}
        S['500  x  300']['8.0'] = RHS(3.000e-01, 5.000e-01, 8.000e-03, 1.200e-02, area=1.250e-02, Iyy=4.370e-04, Izz=2.000e-04, J=4.260e-04)
        S['500  x  300']['8.8'] = RHS(3.000e-01, 5.000e-01, 8.800e-03, 1.320e-02, area=1.370e-02, Iyy=4.780e-04, Izz=2.180e-04, J=4.660e-04)
        S['500  x  300']['10.0'] = RHS(3.000e-01, 5.000e-01, 1.000e-02, 1.500e-02, area=1.550e-02, Iyy=5.380e-04, Izz=2.440e-04, J=5.240e-04)
        S['500  x  300']['11.0'] = RHS(3.000e-01, 5.000e-01, 1.100e-02, 1.650e-02, area=1.700e-02, Iyy=5.860e-04, Izz=2.660e-04, J=5.730e-04)
        S['500  x  300']['12.5'] = RHS(3.000e-01, 5.000e-01, 1.250e-02, 1.875e-02, area=1.920e-02, Iyy=6.580e-04, Izz=2.980e-04, J=6.440e-04)
        S['500  x  300']['14.2'] = RHS(3.000e-01, 5.000e-01, 1.420e-02, 2.130e-02, area=2.170e-02, Iyy=7.370e-04, Izz=3.320e-04, J=7.220e-04)
        S['500  x  300']['16.0'] = RHS(3.000e-01, 5.000e-01, 1.600e-02, 2.400e-02, area=2.430e-02, Iyy=8.180e-04, Izz=3.680e-04, J=8.030e-04)

        return S # No. data rows = 312

    @staticmethod
    def RHS_HotFinished():
        if BlueBook.__RHS_HotFinished is None:
            BlueBook.__RHS_HotFinished = BlueBook.__create_RHS_HotFinished()
        return BlueBook.__RHS_HotFinished

    __SHS_ColdFormed = None

    @staticmethod
    def __create_SHS_ColdFormed():
        S = {}
        S['25  x  25'] = {}
        S['25  x  25']['2.0'] = RHS(2.500e-02, 2.500e-02, 2.000e-03, 5.000e-03, area=1.740e-04, Iyy=1.480e-08, Izz=1.480e-08, J=2.530e-08)
        S['25  x  25']['2.5'] = RHS(2.500e-02, 2.500e-02, 2.500e-03, 6.250e-03, area=2.090e-04, Iyy=1.690e-08, Izz=1.690e-08, J=2.970e-08)
        S['30  x  30'] = {}
        S['30  x  30']['2.0'] = RHS(3.000e-02, 3.000e-02, 2.000e-03, 5.000e-03, area=2.140e-04, Iyy=2.720e-08, Izz=2.720e-08, J=4.540e-08)
        S['30  x  30']['2.5'] = RHS(3.000e-02, 3.000e-02, 2.500e-03, 6.250e-03, area=2.590e-04, Iyy=3.160e-08, Izz=3.160e-08, J=5.400e-08)
        S['30  x  30']['3.0'] = RHS(3.000e-02, 3.000e-02, 3.000e-03, 7.500e-03, area=3.010e-04, Iyy=3.500e-08, Izz=3.500e-08, J=6.150e-08)
        S['40  x  40'] = {}
        S['40  x  40']['2.0'] = RHS(4.000e-02, 4.000e-02, 2.000e-03, 5.000e-03, area=2.940e-04, Iyy=6.940e-08, Izz=6.940e-08, J=1.130e-07)
        S['40  x  40']['2.5'] = RHS(4.000e-02, 4.000e-02, 2.500e-03, 6.250e-03, area=3.590e-04, Iyy=8.220e-08, Izz=8.220e-08, J=1.360e-07)
        S['40  x  40']['3.0'] = RHS(4.000e-02, 4.000e-02, 3.000e-03, 7.500e-03, area=4.210e-04, Iyy=9.320e-08, Izz=9.320e-08, J=1.580e-07)
        S['40  x  40']['4.0'] = RHS(4.000e-02, 4.000e-02, 4.000e-03, 1.000e-02, area=5.350e-04, Iyy=1.110e-07, Izz=1.110e-07, J=1.940e-07)
        S['50  x  50'] = {}
        S['50  x  50']['2.5'] = RHS(5.000e-02, 5.000e-02, 2.500e-03, 6.250e-03, area=4.590e-04, Iyy=1.690e-07, Izz=1.690e-07, J=2.750e-07)
        S['50  x  50']['3.0'] = RHS(5.000e-02, 5.000e-02, 3.000e-03, 7.500e-03, area=5.410e-04, Iyy=1.950e-07, Izz=1.950e-07, J=3.210e-07)
        S['50  x  50']['4.0'] = RHS(5.000e-02, 5.000e-02, 4.000e-03, 1.000e-02, area=6.950e-04, Iyy=2.370e-07, Izz=2.370e-07, J=4.040e-07)
        S['50  x  50']['5.0'] = RHS(5.000e-02, 5.000e-02, 5.000e-03, 1.250e-02, area=8.360e-04, Iyy=2.700e-07, Izz=2.700e-07, J=4.750e-07)
        S['60  x  60'] = {}
        S['60  x  60']['3.0'] = RHS(6.000e-02, 6.000e-02, 3.000e-03, 7.500e-03, area=6.610e-04, Iyy=3.510e-07, Izz=3.510e-07, J=5.710e-07)
        S['60  x  60']['4.0'] = RHS(6.000e-02, 6.000e-02, 4.000e-03, 1.000e-02, area=8.550e-04, Iyy=4.360e-07, Izz=4.360e-07, J=7.260e-07)
        S['60  x  60']['5.0'] = RHS(6.000e-02, 6.000e-02, 5.000e-03, 1.250e-02, area=1.040e-03, Iyy=5.050e-07, Izz=5.050e-07, J=8.640e-07)
        S['60  x  60']['6.0'] = RHS(6.000e-02, 6.000e-02, 6.000e-03, 1.500e-02, area=1.200e-03, Iyy=5.610e-07, Izz=5.610e-07, J=9.840e-07)
        S['70  x  70'] = {}
        S['70  x  70']['3.0'] = RHS(7.000e-02, 7.000e-02, 3.000e-03, 7.500e-03, area=7.810e-04, Iyy=5.750e-07, Izz=5.750e-07, J=9.240e-07)
        S['70  x  70']['3.5'] = RHS(7.000e-02, 7.000e-02, 3.500e-03, 8.750e-03, area=8.990e-04, Iyy=6.510e-07, Izz=6.510e-07, J=1.060e-06)
        S['70  x  70']['4.0'] = RHS(7.000e-02, 7.000e-02, 4.000e-03, 1.000e-02, area=1.010e-03, Iyy=7.210e-07, Izz=7.210e-07, J=1.190e-06)
        S['70  x  70']['5.0'] = RHS(7.000e-02, 7.000e-02, 5.000e-03, 1.250e-02, area=1.240e-03, Iyy=8.460e-07, Izz=8.460e-07, J=1.420e-06)
        S['70  x  70']['6.0'] = RHS(7.000e-02, 7.000e-02, 6.000e-03, 1.500e-02, area=1.440e-03, Iyy=9.520e-07, Izz=9.520e-07, J=1.630e-06)
        S['80  x  80'] = {}
        S['80  x  80']['3.0'] = RHS(8.000e-02, 8.000e-02, 3.000e-03, 7.500e-03, area=9.010e-04, Iyy=8.780e-07, Izz=8.780e-07, J=1.400e-06)
        S['80  x  80']['3.5'] = RHS(8.000e-02, 8.000e-02, 3.500e-03, 8.750e-03, area=1.040e-03, Iyy=9.980e-07, Izz=9.980e-07, J=1.610e-06)
        S['80  x  80']['4.0'] = RHS(8.000e-02, 8.000e-02, 4.000e-03, 1.000e-02, area=1.170e-03, Iyy=1.110e-06, Izz=1.110e-06, J=1.800e-06)
        S['80  x  80']['5.0'] = RHS(8.000e-02, 8.000e-02, 5.000e-03, 1.250e-02, area=1.440e-03, Iyy=1.310e-06, Izz=1.310e-06, J=2.180e-06)
        S['80  x  80']['6.0'] = RHS(8.000e-02, 8.000e-02, 6.000e-03, 1.500e-02, area=1.680e-03, Iyy=1.490e-06, Izz=1.490e-06, J=2.520e-06)
        S['90  x  90'] = {}
        S['90  x  90']['3.0'] = RHS(9.000e-02, 9.000e-02, 3.000e-03, 7.500e-03, area=1.020e-03, Iyy=1.270e-06, Izz=1.270e-06, J=2.010e-06)
        S['90  x  90']['3.5'] = RHS(9.000e-02, 9.000e-02, 3.500e-03, 8.750e-03, area=1.180e-03, Iyy=1.450e-06, Izz=1.450e-06, J=2.320e-06)
        S['90  x  90']['4.0'] = RHS(9.000e-02, 9.000e-02, 4.000e-03, 1.000e-02, area=1.330e-03, Iyy=1.620e-06, Izz=1.620e-06, J=2.610e-06)
        S['90  x  90']['5.0'] = RHS(9.000e-02, 9.000e-02, 5.000e-03, 1.250e-02, area=1.640e-03, Iyy=1.930e-06, Izz=1.930e-06, J=3.160e-06)
        S['90  x  90']['6.0'] = RHS(9.000e-02, 9.000e-02, 6.000e-03, 1.500e-02, area=1.920e-03, Iyy=2.200e-06, Izz=2.200e-06, J=3.680e-06)
        S['100  x  100'] = {}
        S['100  x  100']['3.0'] = RHS(1.000e-01, 1.000e-01, 3.000e-03, 7.500e-03, area=1.140e-03, Iyy=1.770e-06, Izz=1.770e-06, J=2.790e-06)
        S['100  x  100']['4.0'] = RHS(1.000e-01, 1.000e-01, 4.000e-03, 1.000e-02, area=1.490e-03, Iyy=2.260e-06, Izz=2.260e-06, J=3.620e-06)
        S['100  x  100']['5.0'] = RHS(1.000e-01, 1.000e-01, 5.000e-03, 1.250e-02, area=1.840e-03, Iyy=2.710e-06, Izz=2.710e-06, J=4.410e-06)
        S['100  x  100']['6.0'] = RHS(1.000e-01, 1.000e-01, 6.000e-03, 1.500e-02, area=2.160e-03, Iyy=3.110e-06, Izz=3.110e-06, J=5.140e-06)
        S['100  x  100']['8.0'] = RHS(1.000e-01, 1.000e-01, 8.000e-03, 2.000e-02, area=2.720e-03, Iyy=3.660e-06, Izz=3.660e-06, J=6.450e-06)
        S['120  x  120'] = {}
        S['120  x  120']['3.0'] = RHS(1.200e-01, 1.200e-01, 3.000e-03, 7.500e-03, area=1.380e-03, Iyy=3.120e-06, Izz=3.120e-06, J=4.880e-06)
        S['120  x  120']['4.0'] = RHS(1.200e-01, 1.200e-01, 4.000e-03, 1.000e-02, area=1.810e-03, Iyy=4.020e-06, Izz=4.020e-06, J=6.370e-06)
        S['120  x  120']['5.0'] = RHS(1.200e-01, 1.200e-01, 5.000e-03, 1.250e-02, area=2.240e-03, Iyy=4.850e-06, Izz=4.850e-06, J=7.780e-06)
        S['120  x  120']['6.0'] = RHS(1.200e-01, 1.200e-01, 6.000e-03, 1.500e-02, area=2.640e-03, Iyy=5.620e-06, Izz=5.620e-06, J=9.130e-06)
        S['120  x  120']['8.0'] = RHS(1.200e-01, 1.200e-01, 8.000e-03, 2.000e-02, area=3.360e-03, Iyy=6.770e-06, Izz=6.770e-06, J=1.160e-05)
        S['120  x  120']['10.0'] = RHS(1.200e-01, 1.200e-01, 1.000e-02, 2.500e-02, area=4.060e-03, Iyy=7.770e-06, Izz=7.770e-06, J=1.380e-05)
        S['140  x  140'] = {}
        S['140  x  140']['4.0'] = RHS(1.400e-01, 1.400e-01, 4.000e-03, 1.000e-02, area=2.130e-03, Iyy=6.520e-06, Izz=6.520e-06, J=1.020e-05)
        S['140  x  140']['5.0'] = RHS(1.400e-01, 1.400e-01, 5.000e-03, 1.250e-02, area=2.640e-03, Iyy=7.910e-06, Izz=7.910e-06, J=1.260e-05)
        S['140  x  140']['6.0'] = RHS(1.400e-01, 1.400e-01, 6.000e-03, 1.500e-02, area=3.120e-03, Iyy=9.200e-06, Izz=9.200e-06, J=1.480e-05)
        S['140  x  140']['8.0'] = RHS(1.400e-01, 1.400e-01, 8.000e-03, 2.000e-02, area=4.000e-03, Iyy=1.130e-05, Izz=1.130e-05, J=1.900e-05)
        S['140  x  140']['10.0'] = RHS(1.400e-01, 1.400e-01, 1.000e-02, 2.500e-02, area=4.860e-03, Iyy=1.310e-05, Izz=1.310e-05, J=2.270e-05)
        S['150  x  150'] = {}
        S['150  x  150']['4.0'] = RHS(1.500e-01, 1.500e-01, 4.000e-03, 1.000e-02, area=2.290e-03, Iyy=8.080e-06, Izz=8.080e-06, J=1.260e-05)
        S['150  x  150']['5.0'] = RHS(1.500e-01, 1.500e-01, 5.000e-03, 1.250e-02, area=2.840e-03, Iyy=9.820e-06, Izz=9.820e-06, J=1.550e-05)
        S['150  x  150']['6.0'] = RHS(1.500e-01, 1.500e-01, 6.000e-03, 1.500e-02, area=3.360e-03, Iyy=1.150e-05, Izz=1.150e-05, J=1.830e-05)
        S['150  x  150']['8.0'] = RHS(1.500e-01, 1.500e-01, 8.000e-03, 2.000e-02, area=4.320e-03, Iyy=1.410e-05, Izz=1.410e-05, J=2.360e-05)
        S['150  x  150']['10.0'] = RHS(1.500e-01, 1.500e-01, 1.000e-02, 2.500e-02, area=5.260e-03, Iyy=1.650e-05, Izz=1.650e-05, J=2.840e-05)
        S['160  x  160'] = {}
        S['160  x  160']['4.0'] = RHS(1.600e-01, 1.600e-01, 4.000e-03, 1.000e-02, area=2.450e-03, Iyy=9.870e-06, Izz=9.870e-06, J=1.540e-05)
        S['160  x  160']['5.0'] = RHS(1.600e-01, 1.600e-01, 5.000e-03, 1.250e-02, area=3.040e-03, Iyy=1.200e-05, Izz=1.200e-05, J=1.900e-05)
        S['160  x  160']['6.0'] = RHS(1.600e-01, 1.600e-01, 6.000e-03, 1.500e-02, area=3.600e-03, Iyy=1.400e-05, Izz=1.400e-05, J=2.240e-05)
        S['160  x  160']['8.0'] = RHS(1.600e-01, 1.600e-01, 8.000e-03, 2.000e-02, area=4.640e-03, Iyy=1.740e-05, Izz=1.740e-05, J=2.900e-05)
        S['160  x  160']['10.0'] = RHS(1.600e-01, 1.600e-01, 1.000e-02, 2.500e-02, area=5.660e-03, Iyy=2.050e-05, Izz=2.050e-05, J=3.490e-05)
        S['180  x  180'] = {}
        S['180  x  180']['5.0'] = RHS(1.800e-01, 1.800e-01, 5.000e-03, 1.250e-02, area=3.440e-03, Iyy=1.740e-05, Izz=1.740e-05, J=2.720e-05)
        S['180  x  180']['6.0'] = RHS(1.800e-01, 1.800e-01, 6.000e-03, 1.500e-02, area=4.080e-03, Iyy=2.040e-05, Izz=2.040e-05, J=3.220e-05)
        S['180  x  180']['6.3'] = RHS(1.800e-01, 1.800e-01, 6.300e-03, 1.575e-02, area=4.240e-03, Iyy=2.100e-05, Izz=2.100e-05, J=3.380e-05)
        S['180  x  180']['8.0'] = RHS(1.800e-01, 1.800e-01, 8.000e-03, 2.000e-02, area=5.280e-03, Iyy=2.550e-05, Izz=2.550e-05, J=4.190e-05)
        S['180  x  180']['10.0'] = RHS(1.800e-01, 1.800e-01, 1.000e-02, 2.500e-02, area=6.460e-03, Iyy=3.020e-05, Izz=3.020e-05, J=5.070e-05)
        S['180  x  180']['12.0'] = RHS(1.800e-01, 1.800e-01, 1.200e-02, 3.000e-02, area=7.450e-03, Iyy=3.320e-05, Izz=3.320e-05, J=5.860e-05)
        S['180  x  180']['12.5'] = RHS(1.800e-01, 1.800e-01, 1.250e-02, 3.125e-02, area=7.700e-03, Iyy=3.410e-05, Izz=3.410e-05, J=6.050e-05)
        S['200  x  200'] = {}
        S['200  x  200']['5.0'] = RHS(2.000e-01, 2.000e-01, 5.000e-03, 1.250e-02, area=3.840e-03, Iyy=2.410e-05, Izz=2.410e-05, J=3.760e-05)
        S['200  x  200']['6.0'] = RHS(2.000e-01, 2.000e-01, 6.000e-03, 1.500e-02, area=4.560e-03, Iyy=2.830e-05, Izz=2.830e-05, J=4.460e-05)
        S['200  x  200']['6.3'] = RHS(2.000e-01, 2.000e-01, 6.300e-03, 1.575e-02, area=4.740e-03, Iyy=2.920e-05, Izz=2.920e-05, J=4.680e-05)
        S['200  x  200']['8.0'] = RHS(2.000e-01, 2.000e-01, 8.000e-03, 2.000e-02, area=5.920e-03, Iyy=3.570e-05, Izz=3.570e-05, J=5.820e-05)
        S['200  x  200']['10.0'] = RHS(2.000e-01, 2.000e-01, 1.000e-02, 2.500e-02, area=7.260e-03, Iyy=4.250e-05, Izz=4.250e-05, J=7.070e-05)
        S['200  x  200']['12.0'] = RHS(2.000e-01, 2.000e-01, 1.200e-02, 3.000e-02, area=8.410e-03, Iyy=4.730e-05, Izz=4.730e-05, J=8.230e-05)
        S['200  x  200']['12.5'] = RHS(2.000e-01, 2.000e-01, 1.250e-02, 3.125e-02, area=8.700e-03, Iyy=4.860e-05, Izz=4.860e-05, J=8.500e-05)
        S['250  x  250'] = {}
        S['250  x  250']['6.0'] = RHS(2.500e-01, 2.500e-01, 6.000e-03, 1.500e-02, area=5.760e-03, Iyy=5.670e-05, Izz=5.670e-05, J=8.840e-05)
        S['250  x  250']['6.3'] = RHS(2.500e-01, 2.500e-01, 6.300e-03, 1.575e-02, area=6.000e-03, Iyy=5.870e-05, Izz=5.870e-05, J=9.290e-05)
        S['250  x  250']['8.0'] = RHS(2.500e-01, 2.500e-01, 8.000e-03, 2.000e-02, area=7.520e-03, Iyy=7.230e-05, Izz=7.230e-05, J=1.160e-04)
        S['250  x  250']['10.0'] = RHS(2.500e-01, 2.500e-01, 1.000e-02, 2.500e-02, area=9.260e-03, Iyy=8.710e-05, Izz=8.710e-05, J=1.420e-04)
        S['250  x  250']['12.0'] = RHS(2.500e-01, 2.500e-01, 1.200e-02, 3.000e-02, area=1.080e-02, Iyy=9.860e-05, Izz=9.860e-05, J=1.670e-04)
        S['250  x  250']['12.5'] = RHS(2.500e-01, 2.500e-01, 1.250e-02, 3.125e-02, area=1.120e-02, Iyy=1.020e-04, Izz=1.020e-04, J=1.730e-04)
        S['300  x  300'] = {}
        S['300  x  300']['6.0'] = RHS(3.000e-01, 3.000e-01, 6.000e-03, 1.500e-02, area=6.960e-03, Iyy=9.960e-05, Izz=9.960e-05, J=1.540e-04)
        S['300  x  300']['6.3'] = RHS(3.000e-01, 3.000e-01, 6.300e-03, 1.575e-02, area=7.260e-03, Iyy=1.030e-04, Izz=1.030e-04, J=1.620e-04)
        S['300  x  300']['8.0'] = RHS(3.000e-01, 3.000e-01, 8.000e-03, 2.000e-02, area=9.120e-03, Iyy=1.280e-04, Izz=1.280e-04, J=2.030e-04)
        S['300  x  300']['10.0'] = RHS(3.000e-01, 3.000e-01, 1.000e-02, 2.500e-02, area=1.130e-02, Iyy=1.550e-04, Izz=1.550e-04, J=2.500e-04)
        S['300  x  300']['12.0'] = RHS(3.000e-01, 3.000e-01, 1.200e-02, 3.000e-02, area=1.320e-02, Iyy=1.780e-04, Izz=1.780e-04, J=2.950e-04)
        S['300  x  300']['12.5'] = RHS(3.000e-01, 3.000e-01, 1.250e-02, 3.125e-02, area=1.370e-02, Iyy=1.830e-04, Izz=1.830e-04, J=3.060e-04)
        S['350  x  350'] = {}
        S['350  x  350']['6.0'] = RHS(3.500e-01, 3.500e-01, 6.000e-03, 1.500e-02, area=8.160e-03, Iyy=1.600e-04, Izz=1.600e-04, J=2.470e-04)
        S['350  x  350']['6.3'] = RHS(3.500e-01, 3.500e-01, 6.300e-03, 1.575e-02, area=8.520e-03, Iyy=1.660e-04, Izz=1.660e-04, J=2.590e-04)
        S['350  x  350']['8.0'] = RHS(3.500e-01, 3.500e-01, 8.000e-03, 2.000e-02, area=1.070e-02, Iyy=2.070e-04, Izz=2.070e-04, J=3.260e-04)
        S['350  x  350']['10.0'] = RHS(3.500e-01, 3.500e-01, 1.000e-02, 2.500e-02, area=1.330e-02, Iyy=2.520e-04, Izz=2.520e-04, J=4.010e-04)
        S['350  x  350']['12.0'] = RHS(3.500e-01, 3.500e-01, 1.200e-02, 3.000e-02, area=1.560e-02, Iyy=2.910e-04, Izz=2.910e-04, J=4.760e-04)
        S['350  x  350']['12.5'] = RHS(3.500e-01, 3.500e-01, 1.250e-02, 3.125e-02, area=1.620e-02, Iyy=3.000e-04, Izz=3.000e-04, J=4.940e-04)
        S['400  x  400'] = {}
        S['400  x  400']['6.0'] = RHS(4.000e-01, 4.000e-01, 6.000e-03, 1.500e-02, area=9.360e-03, Iyy=2.410e-04, Izz=2.410e-04, J=3.700e-04)
        S['400  x  400']['6.3'] = RHS(4.000e-01, 4.000e-01, 6.300e-03, 1.575e-02, area=9.780e-03, Iyy=2.510e-04, Izz=2.510e-04, J=3.890e-04)
        S['400  x  400']['8.0'] = RHS(4.000e-01, 4.000e-01, 8.000e-03, 2.000e-02, area=1.230e-02, Iyy=3.130e-04, Izz=3.130e-04, J=4.890e-04)
        S['400  x  400']['10.0'] = RHS(4.000e-01, 4.000e-01, 1.000e-02, 2.500e-02, area=1.530e-02, Iyy=3.820e-04, Izz=3.820e-04, J=6.040e-04)
        S['400  x  400']['12.0'] = RHS(4.000e-01, 4.000e-01, 1.200e-02, 3.000e-02, area=1.800e-02, Iyy=4.430e-04, Izz=4.430e-04, J=7.180e-04)
        S['400  x  400']['12.5'] = RHS(4.000e-01, 4.000e-01, 1.250e-02, 3.125e-02, area=1.870e-02, Iyy=4.590e-04, Izz=4.590e-04, J=7.460e-04)

        return S # No. data rows = 96

    @staticmethod
    def SHS_ColdFormed():
        if BlueBook.__SHS_ColdFormed is None:
            BlueBook.__SHS_ColdFormed = BlueBook.__create_SHS_ColdFormed()
        return BlueBook.__SHS_ColdFormed

    __SHS_HotFinished = None

    @staticmethod
    def __create_SHS_HotFinished():
        S = {}
        S['40  x  40'] = {}
        S['40  x  40']['2.9'] = RHS(4.000e-02, 4.000e-02, 2.900e-03, 4.350e-03, area=4.210e-04, Iyy=9.540e-08, Izz=9.540e-08, J=1.530e-07)
        S['40  x  40']['3.0'] = RHS(4.000e-02, 4.000e-02, 3.000e-03, 4.500e-03, area=4.340e-04, Iyy=9.780e-08, Izz=9.780e-08, J=1.570e-07)
        S['40  x  40']['3.2'] = RHS(4.000e-02, 4.000e-02, 3.200e-03, 4.800e-03, area=4.600e-04, Iyy=1.020e-07, Izz=1.020e-07, J=1.650e-07)
        S['40  x  40']['3.6'] = RHS(4.000e-02, 4.000e-02, 3.600e-03, 5.400e-03, area=5.100e-04, Iyy=1.110e-07, Izz=1.110e-07, J=1.810e-07)
        S['40  x  40']['4.0'] = RHS(4.000e-02, 4.000e-02, 4.000e-03, 6.000e-03, area=5.590e-04, Iyy=1.180e-07, Izz=1.180e-07, J=1.950e-07)
        S['40  x  40']['5.0'] = RHS(4.000e-02, 4.000e-02, 5.000e-03, 7.500e-03, area=6.730e-04, Iyy=1.340e-07, Izz=1.340e-07, J=2.250e-07)
        S['40  x  40']['5.6'] = RHS(4.000e-02, 4.000e-02, 5.600e-03, 8.400e-03, area=7.370e-04, Iyy=1.410e-07, Izz=1.410e-07, J=2.400e-07)
        S['40  x  40']['6.3'] = RHS(4.000e-02, 4.000e-02, 6.300e-03, 9.450e-03, area=8.070e-04, Iyy=1.470e-07, Izz=1.470e-07, J=2.540e-07)
        S['50  x  50'] = {}
        S['50  x  50']['2.9'] = RHS(5.000e-02, 5.000e-02, 2.900e-03, 4.350e-03, area=5.370e-04, Iyy=1.970e-07, Izz=1.970e-07, J=3.120e-07)
        S['50  x  50']['3.0'] = RHS(5.000e-02, 5.000e-02, 3.000e-03, 4.500e-03, area=5.540e-04, Iyy=2.020e-07, Izz=2.020e-07, J=3.210e-07)
        S['50  x  50']['3.2'] = RHS(5.000e-02, 5.000e-02, 3.200e-03, 4.800e-03, area=5.880e-04, Iyy=2.120e-07, Izz=2.120e-07, J=3.380e-07)
        S['50  x  50']['3.6'] = RHS(5.000e-02, 5.000e-02, 3.600e-03, 5.400e-03, area=6.540e-04, Iyy=2.320e-07, Izz=2.320e-07, J=3.720e-07)
        S['50  x  50']['4.0'] = RHS(5.000e-02, 5.000e-02, 4.000e-03, 6.000e-03, area=7.190e-04, Iyy=2.500e-07, Izz=2.500e-07, J=4.040e-07)
        S['50  x  50']['5.0'] = RHS(5.000e-02, 5.000e-02, 5.000e-03, 7.500e-03, area=8.730e-04, Iyy=2.890e-07, Izz=2.890e-07, J=4.760e-07)
        S['50  x  50']['5.6'] = RHS(5.000e-02, 5.000e-02, 5.600e-03, 8.400e-03, area=9.610e-04, Iyy=3.080e-07, Izz=3.080e-07, J=5.130e-07)
        S['50  x  50']['6.3'] = RHS(5.000e-02, 5.000e-02, 6.300e-03, 9.450e-03, area=1.060e-03, Iyy=3.280e-07, Izz=3.280e-07, J=5.520e-07)
        S['50  x  50']['7.1'] = RHS(5.000e-02, 5.000e-02, 7.100e-03, 1.065e-02, area=1.160e-03, Iyy=3.450e-07, Izz=3.450e-07, J=5.890e-07)
        S['50  x  50']['8.0'] = RHS(5.000e-02, 5.000e-02, 8.000e-03, 1.200e-02, area=1.280e-03, Iyy=3.600e-07, Izz=3.600e-07, J=6.230e-07)
        S['60  x  60'] = {}
        S['60  x  60']['2.9'] = RHS(6.000e-02, 6.000e-02, 2.900e-03, 4.350e-03, area=6.530e-04, Iyy=3.520e-07, Izz=3.520e-07, J=5.530e-07)
        S['60  x  60']['3.0'] = RHS(6.000e-02, 6.000e-02, 3.000e-03, 4.500e-03, area=6.740e-04, Iyy=3.620e-07, Izz=3.620e-07, J=5.690e-07)
        S['60  x  60']['3.2'] = RHS(6.000e-02, 6.000e-02, 3.200e-03, 4.800e-03, area=7.160e-04, Iyy=3.820e-07, Izz=3.820e-07, J=6.020e-07)
        S['60  x  60']['3.6'] = RHS(6.000e-02, 6.000e-02, 3.600e-03, 5.400e-03, area=7.980e-04, Iyy=4.190e-07, Izz=4.190e-07, J=6.650e-07)
        S['60  x  60']['4.0'] = RHS(6.000e-02, 6.000e-02, 4.000e-03, 6.000e-03, area=8.790e-04, Iyy=4.540e-07, Izz=4.540e-07, J=7.250e-07)
        S['60  x  60']['5.0'] = RHS(6.000e-02, 6.000e-02, 5.000e-03, 7.500e-03, area=1.070e-03, Iyy=5.330e-07, Izz=5.330e-07, J=8.640e-07)
        S['60  x  60']['5.6'] = RHS(6.000e-02, 6.000e-02, 5.600e-03, 8.400e-03, area=1.180e-03, Iyy=5.740e-07, Izz=5.740e-07, J=9.390e-07)
        S['60  x  60']['6.3'] = RHS(6.000e-02, 6.000e-02, 6.300e-03, 9.450e-03, area=1.310e-03, Iyy=6.160e-07, Izz=6.160e-07, J=1.020e-06)
        S['60  x  60']['7.1'] = RHS(6.000e-02, 6.000e-02, 7.100e-03, 1.065e-02, area=1.450e-03, Iyy=6.580e-07, Izz=6.580e-07, J=1.100e-06)
        S['60  x  60']['8.0'] = RHS(6.000e-02, 6.000e-02, 8.000e-03, 1.200e-02, area=1.600e-03, Iyy=6.970e-07, Izz=6.970e-07, J=1.180e-06)
        S['70  x  70'] = {}
        S['70  x  70']['3.0'] = RHS(7.000e-02, 7.000e-02, 3.000e-03, 4.500e-03, area=7.940e-04, Iyy=5.900e-07, Izz=5.900e-07, J=9.220e-07)
        S['70  x  70']['3.2'] = RHS(7.000e-02, 7.000e-02, 3.200e-03, 4.800e-03, area=8.440e-04, Iyy=6.230e-07, Izz=6.230e-07, J=9.760e-07)
        S['70  x  70']['3.6'] = RHS(7.000e-02, 7.000e-02, 3.600e-03, 5.400e-03, area=9.420e-04, Iyy=6.860e-07, Izz=6.860e-07, J=1.080e-06)
        S['70  x  70']['4.0'] = RHS(7.000e-02, 7.000e-02, 4.000e-03, 6.000e-03, area=1.040e-03, Iyy=7.470e-07, Izz=7.470e-07, J=1.180e-06)
        S['70  x  70']['5.0'] = RHS(7.000e-02, 7.000e-02, 5.000e-03, 7.500e-03, area=1.270e-03, Iyy=8.850e-07, Izz=8.850e-07, J=1.420e-06)
        S['70  x  70']['5.6'] = RHS(7.000e-02, 7.000e-02, 5.600e-03, 8.400e-03, area=1.410e-03, Iyy=9.590e-07, Izz=9.590e-07, J=1.550e-06)
        S['70  x  70']['6.3'] = RHS(7.000e-02, 7.000e-02, 6.300e-03, 9.450e-03, area=1.560e-03, Iyy=1.040e-06, Izz=1.040e-06, J=1.690e-06)
        S['70  x  70']['7.1'] = RHS(7.000e-02, 7.000e-02, 7.100e-03, 1.065e-02, area=1.730e-03, Iyy=1.120e-06, Izz=1.120e-06, J=1.850e-06)
        S['70  x  70']['8.0'] = RHS(7.000e-02, 7.000e-02, 8.000e-03, 1.200e-02, area=1.920e-03, Iyy=1.200e-06, Izz=1.200e-06, J=2.000e-06)
        S['70  x  70']['8.8'] = RHS(7.000e-02, 7.000e-02, 8.800e-03, 1.320e-02, area=2.070e-03, Iyy=1.260e-06, Izz=1.260e-06, J=2.120e-06)
        S['80  x  80'] = {}
        S['80  x  80']['3.2'] = RHS(8.000e-02, 8.000e-02, 3.200e-03, 4.800e-03, area=9.720e-04, Iyy=9.500e-07, Izz=9.500e-07, J=1.480e-06)
        S['80  x  80']['3.6'] = RHS(8.000e-02, 8.000e-02, 3.600e-03, 5.400e-03, area=1.090e-03, Iyy=1.050e-06, Izz=1.050e-06, J=1.640e-06)
        S['80  x  80']['4.0'] = RHS(8.000e-02, 8.000e-02, 4.000e-03, 6.000e-03, area=1.200e-03, Iyy=1.140e-06, Izz=1.140e-06, J=1.800e-06)
        S['80  x  80']['5.0'] = RHS(8.000e-02, 8.000e-02, 5.000e-03, 7.500e-03, area=1.470e-03, Iyy=1.370e-06, Izz=1.370e-06, J=2.170e-06)
        S['80  x  80']['5.6'] = RHS(8.000e-02, 8.000e-02, 5.600e-03, 8.400e-03, area=1.630e-03, Iyy=1.490e-06, Izz=1.490e-06, J=2.380e-06)
        S['80  x  80']['6.3'] = RHS(8.000e-02, 8.000e-02, 6.300e-03, 9.450e-03, area=1.810e-03, Iyy=1.620e-06, Izz=1.620e-06, J=2.620e-06)
        S['80  x  80']['7.1'] = RHS(8.000e-02, 8.000e-02, 7.100e-03, 1.065e-02, area=2.020e-03, Iyy=1.760e-06, Izz=1.760e-06, J=2.860e-06)
        S['80  x  80']['8.0'] = RHS(8.000e-02, 8.000e-02, 8.000e-03, 1.200e-02, area=2.240e-03, Iyy=1.890e-06, Izz=1.890e-06, J=3.120e-06)
        S['80  x  80']['8.8'] = RHS(8.000e-02, 8.000e-02, 8.800e-03, 1.320e-02, area=2.420e-03, Iyy=2.000e-06, Izz=2.000e-06, J=3.320e-06)
        S['80  x  80']['10.0'] = RHS(8.000e-02, 8.000e-02, 1.000e-02, 1.500e-02, area=2.690e-03, Iyy=2.140e-06, Izz=2.140e-06, J=3.600e-06)
        S['90  x  90'] = {}
        S['90  x  90']['3.6'] = RHS(9.000e-02, 9.000e-02, 3.600e-03, 5.400e-03, area=1.230e-03, Iyy=1.520e-06, Izz=1.520e-06, J=2.370e-06)
        S['90  x  90']['4.0'] = RHS(9.000e-02, 9.000e-02, 4.000e-03, 6.000e-03, area=1.360e-03, Iyy=1.660e-06, Izz=1.660e-06, J=2.600e-06)
        S['90  x  90']['5.0'] = RHS(9.000e-02, 9.000e-02, 5.000e-03, 7.500e-03, area=1.670e-03, Iyy=2.000e-06, Izz=2.000e-06, J=3.160e-06)
        S['90  x  90']['5.6'] = RHS(9.000e-02, 9.000e-02, 5.600e-03, 8.400e-03, area=1.860e-03, Iyy=2.180e-06, Izz=2.180e-06, J=3.470e-06)
        S['90  x  90']['6.3'] = RHS(9.000e-02, 9.000e-02, 6.300e-03, 9.450e-03, area=2.070e-03, Iyy=2.380e-06, Izz=2.380e-06, J=3.820e-06)
        S['90  x  90']['7.1'] = RHS(9.000e-02, 9.000e-02, 7.100e-03, 1.065e-02, area=2.300e-03, Iyy=2.600e-06, Izz=2.600e-06, J=4.190e-06)
        S['90  x  90']['8.0'] = RHS(9.000e-02, 9.000e-02, 8.000e-03, 1.200e-02, area=2.560e-03, Iyy=2.810e-06, Izz=2.810e-06, J=4.590e-06)
        S['90  x  90']['8.8'] = RHS(9.000e-02, 9.000e-02, 8.800e-03, 1.320e-02, area=2.780e-03, Iyy=2.990e-06, Izz=2.990e-06, J=4.920e-06)
        S['90  x  90']['10.0'] = RHS(9.000e-02, 9.000e-02, 1.000e-02, 1.500e-02, area=3.090e-03, Iyy=3.220e-06, Izz=3.220e-06, J=5.360e-06)
        S['100  x  100'] = {}
        S['100  x  100']['3.6'] = RHS(1.000e-01, 1.000e-01, 3.600e-03, 5.400e-03, area=1.370e-03, Iyy=2.120e-06, Izz=2.120e-06, J=3.280e-06)
        S['100  x  100']['4.0'] = RHS(1.000e-01, 1.000e-01, 4.000e-03, 6.000e-03, area=1.520e-03, Iyy=2.320e-06, Izz=2.320e-06, J=3.610e-06)
        S['100  x  100']['5.0'] = RHS(1.000e-01, 1.000e-01, 5.000e-03, 7.500e-03, area=1.870e-03, Iyy=2.790e-06, Izz=2.790e-06, J=4.390e-06)
        S['100  x  100']['5.6'] = RHS(1.000e-01, 1.000e-01, 5.600e-03, 8.400e-03, area=2.080e-03, Iyy=3.060e-06, Izz=3.060e-06, J=4.840e-06)
        S['100  x  100']['6.3'] = RHS(1.000e-01, 1.000e-01, 6.300e-03, 9.450e-03, area=2.320e-03, Iyy=3.360e-06, Izz=3.360e-06, J=5.340e-06)
        S['100  x  100']['7.1'] = RHS(1.000e-01, 1.000e-01, 7.100e-03, 1.065e-02, area=2.580e-03, Iyy=3.670e-06, Izz=3.670e-06, J=5.890e-06)
        S['100  x  100']['8.0'] = RHS(1.000e-01, 1.000e-01, 8.000e-03, 1.200e-02, area=2.880e-03, Iyy=4.000e-06, Izz=4.000e-06, J=6.460e-06)
        S['100  x  100']['8.8'] = RHS(1.000e-01, 1.000e-01, 8.800e-03, 1.320e-02, area=3.130e-03, Iyy=4.260e-06, Izz=4.260e-06, J=6.940e-06)
        S['100  x  100']['10.0'] = RHS(1.000e-01, 1.000e-01, 1.000e-02, 1.500e-02, area=3.490e-03, Iyy=4.620e-06, Izz=4.620e-06, J=7.610e-06)
        S['120  x  120'] = {}
        S['120  x  120']['4.0'] = RHS(1.200e-01, 1.200e-01, 4.000e-03, 6.000e-03, area=1.840e-03, Iyy=4.100e-06, Izz=4.100e-06, J=6.350e-06)
        S['120  x  120']['5.0'] = RHS(1.200e-01, 1.200e-01, 5.000e-03, 7.500e-03, area=2.270e-03, Iyy=4.980e-06, Izz=4.980e-06, J=7.770e-06)
        S['120  x  120']['5.6'] = RHS(1.200e-01, 1.200e-01, 5.600e-03, 8.400e-03, area=2.530e-03, Iyy=5.470e-06, Izz=5.470e-06, J=8.580e-06)
        S['120  x  120']['6.3'] = RHS(1.200e-01, 1.200e-01, 6.300e-03, 9.450e-03, area=2.820e-03, Iyy=6.030e-06, Izz=6.030e-06, J=9.500e-06)
        S['120  x  120']['7.1'] = RHS(1.200e-01, 1.200e-01, 7.100e-03, 1.065e-02, area=3.150e-03, Iyy=6.630e-06, Izz=6.630e-06, J=1.050e-05)
        S['120  x  120']['8.0'] = RHS(1.200e-01, 1.200e-01, 8.000e-03, 1.200e-02, area=3.520e-03, Iyy=7.260e-06, Izz=7.260e-06, J=1.160e-05)
        S['120  x  120']['8.8'] = RHS(1.200e-01, 1.200e-01, 8.800e-03, 1.320e-02, area=3.830e-03, Iyy=7.790e-06, Izz=7.790e-06, J=1.250e-05)
        S['120  x  120']['10.0'] = RHS(1.200e-01, 1.200e-01, 1.000e-02, 1.500e-02, area=4.290e-03, Iyy=8.520e-06, Izz=8.520e-06, J=1.380e-05)
        S['120  x  120']['11.0'] = RHS(1.200e-01, 1.200e-01, 1.100e-02, 1.650e-02, area=4.670e-03, Iyy=9.080e-06, Izz=9.080e-06, J=1.480e-05)
        S['120  x  120']['12.5'] = RHS(1.200e-01, 1.200e-01, 1.250e-02, 1.875e-02, area=5.210e-03, Iyy=9.820e-06, Izz=9.820e-06, J=1.620e-05)
        S['140  x  140'] = {}
        S['140  x  140']['5.0'] = RHS(1.400e-01, 1.400e-01, 5.000e-03, 7.500e-03, area=2.670e-03, Iyy=8.070e-06, Izz=8.070e-06, J=1.250e-05)
        S['140  x  140']['5.6'] = RHS(1.400e-01, 1.400e-01, 5.600e-03, 8.400e-03, area=2.980e-03, Iyy=8.910e-06, Izz=8.910e-06, J=1.390e-05)
        S['140  x  140']['6.3'] = RHS(1.400e-01, 1.400e-01, 6.300e-03, 9.450e-03, area=3.330e-03, Iyy=9.840e-06, Izz=9.840e-06, J=1.540e-05)
        S['140  x  140']['7.1'] = RHS(1.400e-01, 1.400e-01, 7.100e-03, 1.065e-02, area=3.720e-03, Iyy=1.090e-05, Izz=1.090e-05, J=1.710e-05)
        S['140  x  140']['8.0'] = RHS(1.400e-01, 1.400e-01, 8.000e-03, 1.200e-02, area=4.160e-03, Iyy=1.200e-05, Izz=1.200e-05, J=1.890e-05)
        S['140  x  140']['8.8'] = RHS(1.400e-01, 1.400e-01, 8.800e-03, 1.320e-02, area=4.540e-03, Iyy=1.290e-05, Izz=1.290e-05, J=2.050e-05)
        S['140  x  140']['10.0'] = RHS(1.400e-01, 1.400e-01, 1.000e-02, 1.500e-02, area=5.090e-03, Iyy=1.420e-05, Izz=1.420e-05, J=2.270e-05)
        S['140  x  140']['11.0'] = RHS(1.400e-01, 1.400e-01, 1.100e-02, 1.650e-02, area=5.550e-03, Iyy=1.520e-05, Izz=1.520e-05, J=2.450e-05)
        S['140  x  140']['12.5'] = RHS(1.400e-01, 1.400e-01, 1.250e-02, 1.875e-02, area=6.210e-03, Iyy=1.650e-05, Izz=1.650e-05, J=2.700e-05)
        S['140  x  140']['14.2'] = RHS(1.400e-01, 1.400e-01, 1.420e-02, 2.130e-02, area=6.930e-03, Iyy=1.790e-05, Izz=1.790e-05, J=2.950e-05)
        S['140  x  140']['16.0'] = RHS(1.400e-01, 1.400e-01, 1.600e-02, 2.400e-02, area=7.660e-03, Iyy=1.920e-05, Izz=1.920e-05, J=3.200e-05)
        S['150  x  150'] = {}
        S['150  x  150']['5.0'] = RHS(1.500e-01, 1.500e-01, 5.000e-03, 7.500e-03, area=2.870e-03, Iyy=1.000e-05, Izz=1.000e-05, J=1.550e-05)
        S['150  x  150']['5.6'] = RHS(1.500e-01, 1.500e-01, 5.600e-03, 8.400e-03, area=3.200e-03, Iyy=1.110e-05, Izz=1.110e-05, J=1.720e-05)
        S['150  x  150']['6.3'] = RHS(1.500e-01, 1.500e-01, 6.300e-03, 9.450e-03, area=3.580e-03, Iyy=1.220e-05, Izz=1.220e-05, J=1.910e-05)
        S['150  x  150']['7.1'] = RHS(1.500e-01, 1.500e-01, 7.100e-03, 1.065e-02, area=4.000e-03, Iyy=1.350e-05, Izz=1.350e-05, J=2.120e-05)
        S['150  x  150']['8.0'] = RHS(1.500e-01, 1.500e-01, 8.000e-03, 1.200e-02, area=4.480e-03, Iyy=1.490e-05, Izz=1.490e-05, J=2.350e-05)
        S['150  x  150']['8.8'] = RHS(1.500e-01, 1.500e-01, 8.800e-03, 1.320e-02, area=4.890e-03, Iyy=1.610e-05, Izz=1.610e-05, J=2.550e-05)
        S['150  x  150']['10.0'] = RHS(1.500e-01, 1.500e-01, 1.000e-02, 1.500e-02, area=5.490e-03, Iyy=1.770e-05, Izz=1.770e-05, J=2.830e-05)
        S['150  x  150']['11.0'] = RHS(1.500e-01, 1.500e-01, 1.100e-02, 1.650e-02, area=5.990e-03, Iyy=1.900e-05, Izz=1.900e-05, J=3.060e-05)
        S['150  x  150']['12.5'] = RHS(1.500e-01, 1.500e-01, 1.250e-02, 1.875e-02, area=6.710e-03, Iyy=2.080e-05, Izz=2.080e-05, J=3.380e-05)
        S['150  x  150']['14.2'] = RHS(1.500e-01, 1.500e-01, 1.420e-02, 2.130e-02, area=7.500e-03, Iyy=2.260e-05, Izz=2.260e-05, J=3.710e-05)
        S['150  x  150']['16.0'] = RHS(1.500e-01, 1.500e-01, 1.600e-02, 2.400e-02, area=8.300e-03, Iyy=2.430e-05, Izz=2.430e-05, J=4.030e-05)
        S['160  x  160'] = {}
        S['160  x  160']['5.0'] = RHS(1.600e-01, 1.600e-01, 5.000e-03, 7.500e-03, area=3.070e-03, Iyy=1.220e-05, Izz=1.220e-05, J=1.890e-05)
        S['160  x  160']['5.6'] = RHS(1.600e-01, 1.600e-01, 5.600e-03, 8.400e-03, area=3.420e-03, Iyy=1.350e-05, Izz=1.350e-05, J=2.100e-05)
        S['160  x  160']['6.3'] = RHS(1.600e-01, 1.600e-01, 6.300e-03, 9.450e-03, area=3.830e-03, Iyy=1.500e-05, Izz=1.500e-05, J=2.330e-05)
        S['160  x  160']['7.1'] = RHS(1.600e-01, 1.600e-01, 7.100e-03, 1.065e-02, area=4.290e-03, Iyy=1.660e-05, Izz=1.660e-05, J=2.600e-05)
        S['160  x  160']['8.0'] = RHS(1.600e-01, 1.600e-01, 8.000e-03, 1.200e-02, area=4.800e-03, Iyy=1.830e-05, Izz=1.830e-05, J=2.880e-05)
        S['160  x  160']['8.8'] = RHS(1.600e-01, 1.600e-01, 8.800e-03, 1.320e-02, area=5.240e-03, Iyy=1.980e-05, Izz=1.980e-05, J=3.120e-05)
        S['160  x  160']['10.0'] = RHS(1.600e-01, 1.600e-01, 1.000e-02, 1.500e-02, area=5.890e-03, Iyy=2.190e-05, Izz=2.190e-05, J=3.480e-05)
        S['160  x  160']['11.0'] = RHS(1.600e-01, 1.600e-01, 1.100e-02, 1.650e-02, area=6.430e-03, Iyy=2.350e-05, Izz=2.350e-05, J=3.760e-05)
        S['160  x  160']['12.5'] = RHS(1.600e-01, 1.600e-01, 1.250e-02, 1.875e-02, area=7.210e-03, Iyy=2.580e-05, Izz=2.580e-05, J=4.160e-05)
        S['160  x  160']['14.2'] = RHS(1.600e-01, 1.600e-01, 1.420e-02, 2.130e-02, area=8.070e-03, Iyy=2.810e-05, Izz=2.810e-05, J=4.580e-05)
        S['160  x  160']['16.0'] = RHS(1.600e-01, 1.600e-01, 1.600e-02, 2.400e-02, area=8.940e-03, Iyy=3.030e-05, Izz=3.030e-05, J=4.990e-05)
        S['180  x  180'] = {}
        S['180  x  180']['6.3'] = RHS(1.800e-01, 1.800e-01, 6.300e-03, 9.450e-03, area=4.330e-03, Iyy=2.170e-05, Izz=2.170e-05, J=3.360e-05)
        S['180  x  180']['7.1'] = RHS(1.800e-01, 1.800e-01, 7.100e-03, 1.065e-02, area=4.860e-03, Iyy=2.400e-05, Izz=2.400e-05, J=3.740e-05)
        S['180  x  180']['8.0'] = RHS(1.800e-01, 1.800e-01, 8.000e-03, 1.200e-02, area=5.440e-03, Iyy=2.660e-05, Izz=2.660e-05, J=4.160e-05)
        S['180  x  180']['8.8'] = RHS(1.800e-01, 1.800e-01, 8.800e-03, 1.320e-02, area=5.940e-03, Iyy=2.880e-05, Izz=2.880e-05, J=4.520e-05)
        S['180  x  180']['10.0'] = RHS(1.800e-01, 1.800e-01, 1.000e-02, 1.500e-02, area=6.690e-03, Iyy=3.190e-05, Izz=3.190e-05, J=5.050e-05)
        S['180  x  180']['11.0'] = RHS(1.800e-01, 1.800e-01, 1.100e-02, 1.650e-02, area=7.310e-03, Iyy=3.440e-05, Izz=3.440e-05, J=5.470e-05)
        S['180  x  180']['12.5'] = RHS(1.800e-01, 1.800e-01, 1.250e-02, 1.875e-02, area=8.210e-03, Iyy=3.790e-05, Izz=3.790e-05, J=6.070e-05)
        S['180  x  180']['14.2'] = RHS(1.800e-01, 1.800e-01, 1.420e-02, 2.130e-02, area=9.200e-03, Iyy=4.150e-05, Izz=4.150e-05, J=6.710e-05)
        S['180  x  180']['16.0'] = RHS(1.800e-01, 1.800e-01, 1.600e-02, 2.400e-02, area=1.020e-02, Iyy=4.500e-05, Izz=4.500e-05, J=7.340e-05)
        S['200  x  200'] = {}
        S['200  x  200']['5.0'] = RHS(2.000e-01, 2.000e-01, 5.000e-03, 7.500e-03, area=3.870e-03, Iyy=2.440e-05, Izz=2.440e-05, J=3.760e-05)
        S['200  x  200']['5.6'] = RHS(2.000e-01, 2.000e-01, 5.600e-03, 8.400e-03, area=4.320e-03, Iyy=2.710e-05, Izz=2.710e-05, J=4.170e-05)
        S['200  x  200']['6.3'] = RHS(2.000e-01, 2.000e-01, 6.300e-03, 9.450e-03, area=4.840e-03, Iyy=3.010e-05, Izz=3.010e-05, J=4.650e-05)
        S['200  x  200']['7.1'] = RHS(2.000e-01, 2.000e-01, 7.100e-03, 1.065e-02, area=5.420e-03, Iyy=3.340e-05, Izz=3.340e-05, J=5.190e-05)
        S['200  x  200']['8.0'] = RHS(2.000e-01, 2.000e-01, 8.000e-03, 1.200e-02, area=6.080e-03, Iyy=3.710e-05, Izz=3.710e-05, J=5.780e-05)
        S['200  x  200']['8.8'] = RHS(2.000e-01, 2.000e-01, 8.800e-03, 1.320e-02, area=6.650e-03, Iyy=4.020e-05, Izz=4.020e-05, J=6.290e-05)
        S['200  x  200']['10.0'] = RHS(2.000e-01, 2.000e-01, 1.000e-02, 1.500e-02, area=7.490e-03, Iyy=4.470e-05, Izz=4.470e-05, J=7.030e-05)
        S['200  x  200']['11.0'] = RHS(2.000e-01, 2.000e-01, 1.100e-02, 1.650e-02, area=8.190e-03, Iyy=4.830e-05, Izz=4.830e-05, J=7.630e-05)
        S['200  x  200']['12.5'] = RHS(2.000e-01, 2.000e-01, 1.250e-02, 1.875e-02, area=9.210e-03, Iyy=5.340e-05, Izz=5.340e-05, J=8.490e-05)
        S['200  x  200']['14.2'] = RHS(2.000e-01, 2.000e-01, 1.420e-02, 2.130e-02, area=1.030e-02, Iyy=5.870e-05, Izz=5.870e-05, J=9.420e-05)
        S['200  x  200']['16.0'] = RHS(2.000e-01, 2.000e-01, 1.600e-02, 2.400e-02, area=1.150e-02, Iyy=6.390e-05, Izz=6.390e-05, J=1.030e-04)
        S['220  x  220'] = {}
        S['220  x  220']['8.0'] = RHS(2.200e-01, 2.200e-01, 8.000e-03, 1.200e-02, area=6.720e-03, Iyy=5.000e-05, Izz=5.000e-05, J=7.760e-05)
        S['220  x  220']['8.8'] = RHS(2.200e-01, 2.200e-01, 8.800e-03, 1.320e-02, area=7.350e-03, Iyy=5.430e-05, Izz=5.430e-05, J=8.460e-05)
        S['220  x  220']['10.0'] = RHS(2.200e-01, 2.200e-01, 1.000e-02, 1.500e-02, area=8.290e-03, Iyy=6.050e-05, Izz=6.050e-05, J=9.470e-05)
        S['220  x  220']['11.0'] = RHS(2.200e-01, 2.200e-01, 1.100e-02, 1.650e-02, area=9.070e-03, Iyy=6.550e-05, Izz=6.550e-05, J=1.030e-04)
        S['220  x  220']['12.5'] = RHS(2.200e-01, 2.200e-01, 1.250e-02, 1.875e-02, area=1.020e-02, Iyy=7.250e-05, Izz=7.250e-05, J=1.150e-04)
        S['220  x  220']['14.2'] = RHS(2.200e-01, 2.200e-01, 1.420e-02, 2.130e-02, area=1.150e-02, Iyy=8.010e-05, Izz=8.010e-05, J=1.280e-04)
        S['220  x  220']['16.0'] = RHS(2.200e-01, 2.200e-01, 1.600e-02, 2.400e-02, area=1.280e-02, Iyy=8.750e-05, Izz=8.750e-05, J=1.410e-04)
        S['250  x  250'] = {}
        S['250  x  250']['5.0'] = RHS(2.500e-01, 2.500e-01, 5.000e-03, 7.500e-03, area=4.870e-03, Iyy=4.860e-05, Izz=4.860e-05, J=7.430e-05)
        S['250  x  250']['5.6'] = RHS(2.500e-01, 2.500e-01, 5.600e-03, 8.400e-03, area=5.440e-03, Iyy=5.400e-05, Izz=5.400e-05, J=8.270e-05)
        S['250  x  250']['6.3'] = RHS(2.500e-01, 2.500e-01, 6.300e-03, 9.450e-03, area=6.100e-03, Iyy=6.010e-05, Izz=6.010e-05, J=9.240e-05)
        S['250  x  250']['7.1'] = RHS(2.500e-01, 2.500e-01, 7.100e-03, 1.065e-02, area=6.840e-03, Iyy=6.700e-05, Izz=6.700e-05, J=1.030e-04)
        S['250  x  250']['8.0'] = RHS(2.500e-01, 2.500e-01, 8.000e-03, 1.200e-02, area=7.680e-03, Iyy=7.460e-05, Izz=7.460e-05, J=1.150e-04)
        S['250  x  250']['8.8'] = RHS(2.500e-01, 2.500e-01, 8.800e-03, 1.320e-02, area=8.410e-03, Iyy=8.110e-05, Izz=8.110e-05, J=1.260e-04)
        S['250  x  250']['10.0'] = RHS(2.500e-01, 2.500e-01, 1.000e-02, 1.500e-02, area=9.490e-03, Iyy=9.060e-05, Izz=9.060e-05, J=1.410e-04)
        S['250  x  250']['11.0'] = RHS(2.500e-01, 2.500e-01, 1.100e-02, 1.650e-02, area=1.040e-02, Iyy=9.820e-05, Izz=9.820e-05, J=1.540e-04)
        S['250  x  250']['12.5'] = RHS(2.500e-01, 2.500e-01, 1.250e-02, 1.875e-02, area=1.170e-02, Iyy=1.090e-04, Izz=1.090e-04, J=1.720e-04)
        S['250  x  250']['14.2'] = RHS(2.500e-01, 2.500e-01, 1.420e-02, 2.130e-02, area=1.320e-02, Iyy=1.210e-04, Izz=1.210e-04, J=1.910e-04)
        S['250  x  250']['16.0'] = RHS(2.500e-01, 2.500e-01, 1.600e-02, 2.400e-02, area=1.470e-02, Iyy=1.330e-04, Izz=1.330e-04, J=2.110e-04)
        S['260  x  260'] = {}
        S['260  x  260']['7.1'] = RHS(2.600e-01, 2.600e-01, 7.100e-03, 1.065e-02, area=7.130e-03, Iyy=7.570e-05, Izz=7.570e-05, J=1.160e-04)
        S['260  x  260']['8.0'] = RHS(2.600e-01, 2.600e-01, 8.000e-03, 1.200e-02, area=8.000e-03, Iyy=8.420e-05, Izz=8.420e-05, J=1.300e-04)
        S['260  x  260']['8.8'] = RHS(2.600e-01, 2.600e-01, 8.800e-03, 1.320e-02, area=8.760e-03, Iyy=9.160e-05, Izz=9.160e-05, J=1.420e-04)
        S['260  x  260']['10.0'] = RHS(2.600e-01, 2.600e-01, 1.000e-02, 1.500e-02, area=9.890e-03, Iyy=1.020e-04, Izz=1.020e-04, J=1.590e-04)
        S['260  x  260']['11.0'] = RHS(2.600e-01, 2.600e-01, 1.100e-02, 1.650e-02, area=1.080e-02, Iyy=1.110e-04, Izz=1.110e-04, J=1.730e-04)
        S['260  x  260']['12.5'] = RHS(2.600e-01, 2.600e-01, 1.250e-02, 1.875e-02, area=1.220e-02, Iyy=1.240e-04, Izz=1.240e-04, J=1.940e-04)
        S['260  x  260']['14.2'] = RHS(2.600e-01, 2.600e-01, 1.420e-02, 2.130e-02, area=1.370e-02, Iyy=1.370e-04, Izz=1.370e-04, J=2.170e-04)
        S['260  x  260']['16.0'] = RHS(2.600e-01, 2.600e-01, 1.600e-02, 2.400e-02, area=1.530e-02, Iyy=1.510e-04, Izz=1.510e-04, J=2.390e-04)
        S['300  x  300'] = {}
        S['300  x  300']['6.3'] = RHS(3.000e-01, 3.000e-01, 6.300e-03, 9.450e-03, area=7.360e-03, Iyy=1.050e-04, Izz=1.050e-04, J=1.610e-04)
        S['300  x  300']['7.1'] = RHS(3.000e-01, 3.000e-01, 7.100e-03, 1.065e-02, area=8.260e-03, Iyy=1.180e-04, Izz=1.180e-04, J=1.810e-04)
        S['300  x  300']['8.0'] = RHS(3.000e-01, 3.000e-01, 8.000e-03, 1.200e-02, area=9.280e-03, Iyy=1.310e-04, Izz=1.310e-04, J=2.020e-04)
        S['300  x  300']['8.8'] = RHS(3.000e-01, 3.000e-01, 8.800e-03, 1.320e-02, area=1.020e-02, Iyy=1.430e-04, Izz=1.430e-04, J=2.210e-04)
        S['300  x  300']['10.0'] = RHS(3.000e-01, 3.000e-01, 1.000e-02, 1.500e-02, area=1.150e-02, Iyy=1.600e-04, Izz=1.600e-04, J=2.480e-04)
        S['300  x  300']['11.0'] = RHS(3.000e-01, 3.000e-01, 1.100e-02, 1.650e-02, area=1.260e-02, Iyy=1.740e-04, Izz=1.740e-04, J=2.700e-04)
        S['300  x  300']['12.5'] = RHS(3.000e-01, 3.000e-01, 1.250e-02, 1.875e-02, area=1.420e-02, Iyy=1.940e-04, Izz=1.940e-04, J=3.030e-04)
        S['300  x  300']['14.2'] = RHS(3.000e-01, 3.000e-01, 1.420e-02, 2.130e-02, area=1.600e-02, Iyy=2.160e-04, Izz=2.160e-04, J=3.390e-04)
        S['300  x  300']['16.0'] = RHS(3.000e-01, 3.000e-01, 1.600e-02, 2.400e-02, area=1.790e-02, Iyy=2.380e-04, Izz=2.380e-04, J=3.760e-04)
        S['350  x  350'] = {}
        S['350  x  350']['8.0'] = RHS(3.500e-01, 3.500e-01, 8.000e-03, 1.200e-02, area=1.090e-02, Iyy=2.110e-04, Izz=2.110e-04, J=3.240e-04)
        S['350  x  350']['8.8'] = RHS(3.500e-01, 3.500e-01, 8.800e-03, 1.320e-02, area=1.190e-02, Iyy=2.310e-04, Izz=2.310e-04, J=3.540e-04)
        S['350  x  350']['10.0'] = RHS(3.500e-01, 3.500e-01, 1.000e-02, 1.500e-02, area=1.350e-02, Iyy=2.590e-04, Izz=2.590e-04, J=3.990e-04)
        S['350  x  350']['11.0'] = RHS(3.500e-01, 3.500e-01, 1.100e-02, 1.650e-02, area=1.480e-02, Iyy=2.820e-04, Izz=2.820e-04, J=4.350e-04)
        S['350  x  350']['12.5'] = RHS(3.500e-01, 3.500e-01, 1.250e-02, 1.875e-02, area=1.670e-02, Iyy=3.150e-04, Izz=3.150e-04, J=4.890e-04)
        S['350  x  350']['14.2'] = RHS(3.500e-01, 3.500e-01, 1.420e-02, 2.130e-02, area=1.890e-02, Iyy=3.520e-04, Izz=3.520e-04, J=5.490e-04)
        S['350  x  350']['16.0'] = RHS(3.500e-01, 3.500e-01, 1.600e-02, 2.400e-02, area=2.110e-02, Iyy=3.890e-04, Izz=3.890e-04, J=6.100e-04)
        S['400  x  400'] = {}
        S['400  x  400']['8.0'] = RHS(4.000e-01, 4.000e-01, 8.000e-03, 1.200e-02, area=1.250e-02, Iyy=3.190e-04, Izz=3.190e-04, J=4.870e-04)
        S['400  x  400']['8.8'] = RHS(4.000e-01, 4.000e-01, 8.800e-03, 1.320e-02, area=1.370e-02, Iyy=3.480e-04, Izz=3.480e-04, J=5.330e-04)
        S['400  x  400']['10.0'] = RHS(4.000e-01, 4.000e-01, 1.000e-02, 1.500e-02, area=1.550e-02, Iyy=3.910e-04, Izz=3.910e-04, J=6.010e-04)
        S['400  x  400']['11.0'] = RHS(4.000e-01, 4.000e-01, 1.100e-02, 1.650e-02, area=1.700e-02, Iyy=4.270e-04, Izz=4.270e-04, J=6.570e-04)
        S['400  x  400']['12.5'] = RHS(4.000e-01, 4.000e-01, 1.250e-02, 1.875e-02, area=1.920e-02, Iyy=4.780e-04, Izz=4.780e-04, J=7.390e-04)
        S['400  x  400']['14.2'] = RHS(4.000e-01, 4.000e-01, 1.420e-02, 2.130e-02, area=2.170e-02, Iyy=5.350e-04, Izz=5.350e-04, J=8.300e-04)
        S['400  x  400']['16.0'] = RHS(4.000e-01, 4.000e-01, 1.600e-02, 2.400e-02, area=2.430e-02, Iyy=5.930e-04, Izz=5.930e-04, J=9.240e-04)

        return S # No. data rows = 178

    @staticmethod
    def SHS_HotFinished():
        if BlueBook.__SHS_HotFinished is None:
            BlueBook.__SHS_HotFinished = BlueBook.__create_SHS_HotFinished()
        return BlueBook.__SHS_HotFinished
