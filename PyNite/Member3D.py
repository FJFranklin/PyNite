# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 18:04:56 2017

@author: D. Craig Brinck, SE
"""
# %%
import numpy as np
from numpy import zeros, matrix, transpose, add, subtract, matmul, insert
from numpy.linalg import inv
from PyNite.BeamSegment import BeamSegment
import PyNite.FixedEndReactions

# %%
class MemberBasis(object):
    """
    3D Vector Basis for a beam, with origin at start of beam, i-axis along the beam.
    An optional reference point in the ik-plane (non-colinear with the beam) can be
    specified to orientate the j- & k-axes of the beam.

    Attributes
    ----------
    iNode : Node
        start of the beam, i.e., the origin of the local coordinate system
    jNode : Node
        end of beam, which defines basis vector i
    L : number
        length of beam
    origin : 3D coordinate (numpy 1D array)
        the origin of the local coordinate system
    e_i : 3D vector (numpy 1D array)
        basis vector i
    e_j : 3D vector (numpy 1D array)
        basis vector j
    e_k : 3D vector (numpy 1D array)
        basis vector k
    T : 12x12 matrix (numpy 2D array)
        transformation matrix built from direction cosines
    T_inverse : 12x12 matrix (numpy 2D array)
        inverse (transpose) of the transformation matrix
    """
#%%
    @staticmethod
    def NodeToCoordinate(node):
        """
        Convert PyNite Node to 3D coordinate (numpy 1D array)

        Parameters
        ----------
        node : PyNite Node
            coordinate in space

        Returns
        -------
        coord : 3D coordinate (numpy 1D array)
            coordinate in space
        """
        return np.asarray([node.X, node.Y, node.Z], dtype=np.float64)

#%%
    def __init__(self, iNode, jNode, ik_ref=None): # ik_ref must not be colinear with the beam
        """
        Parameters
        ----------
        iNode : Node
            start of the beam, i.e., the origin of the local coordinate system
        jNode : Node
            end of beam, which defines basis vector i
        ik_ref : 3D coordinate (numpy 1D array), optional
            reference point in ik-plane of the beam (default is None)
        """
        self.iNode = iNode
        self.jNode = jNode

        self.origin = MemberBasis.NodeToCoordinate(iNode)

        self.e_i  = MemberBasis.NodeToCoordinate(jNode) - self.origin
        self.L    = np.linalg.norm(self.e_i)  # beam length
        self.e_i /= self.L

        if ik_ref is None:              # this replicates PyNite default
            if 1 - abs(self.e_i[2]) < 1E-8: # vertical beam, or almost
                ref = np.asarray([self.e_i[2],0,0], dtype=np.float64)
            else:
                ref = np.asarray([0,0,-1], dtype=np.float64)
        else:
            ref = np.asarray(ik_ref, dtype=np.float64) - self.origin

        self.e_j = np.cross(ref, ref - self.e_i)
        self.e_k = np.cross(self.e_i, self.e_j)

        self.e_j /= np.linalg.norm(self.e_j)
        self.e_k /= np.linalg.norm(self.e_k)

        dirCos = [self.e_i, self.e_j, self.e_k]

        self.T = np.zeros((12, 12), dtype=np.float64)
        self.T[0:3,  0:3 ] = dirCos
        self.T[3:6,  3:6 ] = dirCos
        self.T[6:9,  6:9 ] = dirCos
        self.T[9:12, 9:12] = dirCos

        # T is orthonormal, so just transpose to get the inverse
        self.T_inverse = self.T.transpose()

        # Equivalent 4x4 transformation matrix
        #self.T_gl = np.asarray([[*self.e_i,0],[*self.e_j,0],[*self.e_k,0],[*self.origin,1]]).transpose()

#%%
    def ToGlobal(self, coordinate):
        """
        Translate from local coordinate system to global coordinates

        Parameters
        ----------
        coordinate : 3D coordinate (numpy 1D array)
            coordinate in local space

        Returns
        -------
        coord : 3D coordinate (numpy 1D array)
            coordinate in global space
        """
        return self.origin + matmul(coordinate, self.T[0:3,0:3])

# %%
class Member3D(MemberBasis):
    """
    A class representing a 3D frame element in a finite element model.
    """
#%%
    __plt = None # Don't import PyPlot unless and until actually needed

    @staticmethod
    def __axis():
        if Member3D.__plt is None:
            import matplotlib.pyplot as plt
            Member3D.__plt = plt

        fig, ax = Member3D.__plt.subplots()
        return ax

#%%
    def __init__(self, Name, iNode, jNode, material, section, ik_ref=None):
        """
        Initializes a new member.
        """
        MemberBasis.__init__(self, iNode, jNode, ik_ref)
        
        self.Name = Name    # A unique name for the member given by the user
        self.ID = None      # Unique index number for the member assigned by the program
        self.M = material   # Material class object
        self.S = section    # Section class object
        self.PtLoads = []   # A list of point loads & moments applied to the element (Direction, P, x) or (Direction, M, x)
        self.DistLoads = [] # A list of linear distributed loads applied to the element (Direction, w1, w2, x1, x2)
        self.SegmentsZ = [] # A list of mathematically continuous beam segments for z-bending
        self.SegmentsY = [] # A list of mathematically continuous beam segments for y-bending
        self.FER = zeros((12,1)) # The fixed end reaction vector
        self.Releases = [False, False, False, False, False, False, False, False, False, False, False, False]

#%%  
    def Display(self, view, wireframe):
        """
        Displays the member in 3D

        Parameters
        ----------
        view : Viewer3D
            3D viewer for plotting members
        wireframe : boolean
            If true, only plot wireframe
        """
        self.S.Display(view, wireframe, self, self.M.color())

#%%  
    def DisplayResults(self, view, result):
        """
        Displays the member in 3D indicating stress

        Parameters
        ----------
        view : Viewer3D
            3D viewer for plotting members
        result : str
            Stress to display - one of: 'seq' (von Mises equivalent), 'sxx' (axial)
        """
        xcount = 20
        values = np.zeros((xcount,7))

        endvec = self.f()
        Fx = endvec[0,0]
        Mx = endvec[3,0]

        for i in range(0, xcount):
            x = self.L * i / (xcount - 1)

            values[i,0] = x
            values[i,1] = Fx
            values[i,2] = self.Shear('Fy', x)
            values[i,3] = self.Shear('Fz', x)
            values[i,4] = Mx
            values[i,5] = self.Moment('My', x)
            values[i,6] = self.Moment('Mz', x)

        self.S.DisplayResults(view, self, values, result, self.M)

#%%
    def k(self):
        """
        Returns the local condensed stiffness matrix
        """
        
        # Get the local stiffness matrix, partitioned as 4 submatrices in
        # preparation for static condensation
        k11, k12, k21, k22 = self.__k_Partition()
               
        # Calculate the condensed local stiffness matrix
        k_Condensed = subtract(k11, matmul(matmul(k12, inv(k22)), k21))
        
        # Expand the condensed local stiffness matrix
        i=0
        for DOF in self.Releases:
            
            if DOF == True:
                k_Condensed = insert(k_Condensed, i, 0, axis = 0)
                k_Condensed = insert(k_Condensed, i, 0, axis = 1)
                
            i += 1
        
        # Return the local stiffness matrix, with end releases applied
        return k_Condensed
    
#%%
    def __k_Partition(self):
        """
        Partitions the local stiffness matrix in preparation for static
        condensation. Used for applying end releases to the member.
        """
        
        # Get the properties needed to form the local stiffness matrix
        E = self.M.E
        G = self.M.G
        Iy = self.S.Iyy
        Iz = self.S.Izz
        J = self.S.J
        A = self.S.A
        L = self.L
        
        # Create the uncondensed local stiffness matrix
        k = matrix([[A*E/L, 0, 0 , 0, 0, 0, -A*E/L, 0, 0, 0, 0, 0],
                    [0, 12*E*Iz/L**3, 0, 0, 0, 6*E*Iz/L**2, 0, -12*E*Iz/L**3, 0, 0, 0, 6*E*Iz/L**2],
                    [0, 0, 12*E*Iy/L**3, 0, -6*E*Iy/L**2, 0, 0, 0, -12*E*Iy/L**3, 0, -6*E*Iy/L**2, 0],
                    [0, 0, 0, G*J/L, 0, 0, 0, 0, 0, -G*J/L, 0, 0],
                    [0, 0, -6*E*Iy/L**2, 0, 4*E*Iy/L, 0, 0, 0, 6*E*Iy/L**2, 0, 2*E*Iy/L, 0],
                    [0, 6*E*Iz/L**2, 0, 0, 0, 4*E*Iz/L, 0, -6*E*Iz/L**2, 0, 0, 0, 2*E*Iz/L],
                    [-A*E/L, 0, 0, 0, 0, 0, A*E/L, 0, 0, 0, 0, 0],
                    [0, -12*E*Iz/L**3, 0, 0, 0, -6*E*Iz/L**2, 0, 12*E*Iz/L**3, 0, 0, 0, -6*E*Iz/L**2],
                    [0, 0, -12*E*Iy/L**3, 0, 6*E*Iy/L**2, 0, 0, 0, 12*E*Iy/L**3, 0, 6*E*Iy/L**2, 0],
                    [0, 0, 0, -G*J/L, 0, 0, 0, 0, 0, G*J/L, 0, 0],
                    [0, 0, -6*E*Iy/L**2, 0, 2*E*Iy/L, 0, 0, 0, 6*E*Iy/L**2, 0, 4*E*Iy/L, 0],
                    [0, 6*E*Iz/L**2, 0, 0, 0, 2*E*Iz/L, 0, -6*E*Iz/L**2, 0, 0, 0, 4*E*Iz/L]])
        
        # Count the number of released degrees of freedom
        NumReleases = 0
        for DOF in self.Releases:
            if DOF == True:
                NumReleases += 1
                    
        # Initialize each partitioned matrix
        k11 = zeros((12 - NumReleases, 12 - NumReleases))
        k12 = zeros((12 - NumReleases, NumReleases))
        k21 = zeros((NumReleases, 12 - NumReleases))
        k22 = zeros((NumReleases, NumReleases))
        
        # Initialize variables used to track rows and columns as the matrix is partitioned
        m11 = 0
        n11 = 0
        m12 = 0
        n12 = 0
        m21 = 0
        n21 = 0
        m22 = 0
        n22 = 0
        
        # Partition the stiffness matrix
        # Step through each term in the local stiffness matrix (m = row, n = column)
        for m in range(12):
            for n in range(12):
                
                # Determine which partitioned matrix this term belongs in
                if self.Releases[m] == False and self.Releases[n] == False:
                    
                    k11.itemset((m11, n11), k[m, n])
                    
                    n11 += 1
                    if n11 == 12 - NumReleases:
                        n11 = 0
                        m11 += 1            
                        
                elif self.Releases[m] == False and self.Releases[n] == True:
                    
                    k12.itemset((m12, n12), k[m, n])
                    
                    n12 += 1
                    if n12 == NumReleases:
                        n12 = 0
                        m12 += 1
                    
                elif self.Releases[m] == True and self.Releases[n] == False:
                    
                    k21.itemset((m21, n21), k[m, n])
                    
                    n21 += 1
                    if n21 == 12 - NumReleases:
                        n21 = 0
                        m21 += 1
                    
                elif self.Releases[m] == True and self.Releases[n] == True:
                    
                    k22.itemset((m22, n22), k[m, n])
                    
                    n22 += 1
                    if n22 == NumReleases:
                        n22 = 0
                        m22 += 1
        
        # Return the matrix, partitioned into 4 submatrices
        return k11, k12, k21, k22
    
#%%
    def fer(self):
        """
        Returns the member's local fixed end reaction vector
        """
        
        # Partition the local stiffness matrix and local fixed end reaction vector
        k11, k12, k21, k22 = self.__k_Partition()
        fer1, fer2 = self.__fer_Partition()
        
        # Calculate the condensed fixed end reaction vector
        ferCondensed = subtract(fer1, matmul(matmul(k12, inv(k22)), fer2))
        
        # Expand the condensed fixed end reaction vector
        i=0
        for DOF in self.Releases:
            
            if DOF == True:
                ferCondensed = insert(ferCondensed, i, 0, axis = 0)
                
            i += 1
        
        # Return the fixed end reaction vector        
        return ferCondensed
    
#%%
    def __fer_Partition(self):
        """
        Partitions the local stiffness matrix in preparation for static
        condensation. Used for applying end releases to the member.
        """
        
        # Initialize the fixed end reaction vector
        fer = zeros((12,1))
        
        # Sum the fixed end reactions for the point loads & moments
        for ptLoad in self.PtLoads:
            if ptLoad[0] == "Fx":
                fer = add(fer, PyNite.FixedEndReactions.FER_AxialPtLoad(ptLoad[1], ptLoad[2], self.L))
            elif ptLoad[0] == "Fy":
                fer = add(fer, PyNite.FixedEndReactions.FER_PtLoad(ptLoad[1], ptLoad[2], self.L, "Fy"))
            elif ptLoad[0] == "Fz":
                fer = add(fer, PyNite.FixedEndReactions.FER_PtLoad(ptLoad[1], ptLoad[2], self.L, "Fz"))
            elif ptLoad[0] == "Mx":
                fer = fer
            elif ptLoad[0] == "My":
                fer = add(fer, PyNite.FixedEndReactions.FER_Moment(ptLoad[1], ptLoad[2], self.L, "My"))
            elif ptLoad[0] == "Mz":     
                fer = add(fer, PyNite.FixedEndReactions.FER_Moment(ptLoad[1], ptLoad[2], self.L, "Mz"))
                
        # Sum the fixed end reactions for the distributed loads
        for distLoad in self.DistLoads:
            fer = add(fer, PyNite.FixedEndReactions.FER_LinLoad(distLoad[1], distLoad[2], distLoad[3], distLoad[4], self.L, distLoad[0]))        
 
        # Count the number of released degrees of freedom
        NumReleases = 0
        for DOF in self.Releases:
            if DOF == True:
                NumReleases += 1
        
        # Partition the fixed end reaction vector
        # Initialize each partitioned fixed end reaction vector
        fer1 = zeros((12 - NumReleases, 1))
        fer2 = zeros((NumReleases, 1))
        
        # Variables used to track indexing during partitioning
        m1 = 0
        m2 = 0
        
        # Partition the fixed end reaction vector
        for m in range(12):
            
            if self.Releases[m] == False:
                
                fer1.itemset((m1, 0), fer.item(m, 0))
                m1 += 1
            
            elif self.Releases == True:
                
                fer2.itemset((m2, 0), fer.item(m, 0))
                m2 += 1
        
        # Return the vector partitioned as 2 subvectors
        return fer1, fer2
    
#%%
    def __fer_Unc(self):
        """
        Returns the member's local fixed end reaction vector, ignoring the effects of end releases.
        Needed to apply the slope-deflection equation properly.
        """
        
        # Initialize the fixed end reaction vector
        fer = zeros((12,1))
        
        # Sum the fixed end reactions for the point loads & moments
        for ptLoad in self.PtLoads:
            if ptLoad[0] == "Fx":
                fer = add(fer, PyNite.FixedEndReactions.FER_AxialPtLoad(ptLoad[1], ptLoad[2], self.L))
            elif ptLoad[0] == "Fy":
                fer = add(fer, PyNite.FixedEndReactions.FER_PtLoad(ptLoad[1], ptLoad[2], self.L, "Fy"))
            elif ptLoad[0] == "Fz":
                fer = add(fer, PyNite.FixedEndReactions.FER_PtLoad(ptLoad[1], ptLoad[2], self.L, "Fz"))
            elif ptLoad[0] == "Mx":
                fer = fer
            elif ptLoad[0] == "My":
                fer = add(fer, PyNite.FixedEndReactions.FER_Moment(ptLoad[1], ptLoad[2], self.L, "My"))
            elif ptLoad[0] == "Mz":     
                fer = add(fer, PyNite.FixedEndReactions.FER_Moment(ptLoad[1], ptLoad[2], self.L, "Mz"))
                
        # Sum the fixed end reactions for the distributed loads
        for distLoad in self.DistLoads:
            fer = add(fer, PyNite.FixedEndReactions.FER_LinLoad(distLoad[1], distLoad[2], distLoad[3], distLoad[4], self.L, distLoad[0]))
        
        # Return the fixed end reaction vector, uncondensed
        return fer

#%%   
    def f(self):
        """
        Returns the member's local end force vector
        """
        
        # Calculate and return the member's local end force vector
        return add(matmul(self.k(), self.d()), self.fer())

#%%
    def d(self):
       """
       Returns the member's local displacement vector
       """

       # Calculate and return the local displacement vector
       return matmul(self.T, self.D())
        
#%%
    # Member global stiffness matrix
    def K(self):
        
        # Calculate and return the stiffness matrix in global coordinates
        return matmul(matmul(self.T_inverse, self.k()), self.T)

#%%
    def F(self):
        
        # Calculate and return the global force vector
        return matmul(self.T_inverse, self.f())
    
#%% 
    # Global fixed end reaction vector
    def FER(self):
        
        # Calculate and return the fixed end reaction vector
        return matmul(self.T_inverse, self.fer())

#%%
    def D(self):
        """
        Returns the member's global displacement vector
        """
        
        # Initialize the displacement vector
        Dvector = zeros((12,1))
        
        # Read in the global displacements from the nodes
        Dvector.itemset((0, 0), self.iNode.DX)
        Dvector.itemset((1, 0), self.iNode.DY)
        Dvector.itemset((2, 0), self.iNode.DZ)
        Dvector.itemset((3, 0), self.iNode.RX)
        Dvector.itemset((4, 0), self.iNode.RY)
        Dvector.itemset((5, 0), self.iNode.RZ)
        Dvector.itemset((6, 0), self.jNode.DX)
        Dvector.itemset((7, 0), self.jNode.DY)
        Dvector.itemset((8, 0), self.jNode.DZ)
        Dvector.itemset((9, 0), self.jNode.RX)
        Dvector.itemset((10, 0), self.jNode.RY)
        Dvector.itemset((11, 0), self.jNode.RZ)
        
        # Return the global displacement vector
        return Dvector

#%%
    # Adds a concentrated moment to the frame element
    def AddMoment(self, M, x, Direction):
        
        self.Moments.append((M, x, Direction))

#%%
    def Shear(self, Direction, x):
        """
        Returns the shear at a point along the member's length
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the shear. Must be one of the following:
                "Fy" = Shear acting on the local y-axis
                "Fz" = Shear acting on the local z-axis
        x : number
            The location at which to find the shear
        """
        
        # Check which direction is of interest
        if Direction == "Fy":
            
            # Check which segment "x" falls on
            for segment in self.SegmentsZ:
                if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):
                    return segment.Shear(x - segment.x1)
                
            if round(x, 10) == round(self.L, 10):  
                lastIndex = len(self.SegmentsZ) - 1
                return self.SegmentsZ[lastIndex].Shear(x - self.SegmentsZ[lastIndex].x1)
                
        elif Direction == "Fz":
            
            for segment in self.SegmentsY:
                
                if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                    
                    return segment.Shear(x - segment.x1)
                
            if round(x,10) == round(self.L,10):
                
                lastIndex = len(self.SegmentsY) - 1
                return self.SegmentsY[lastIndex].Shear(x - self.SegmentsY[lastIndex].x1)
            
#%%
    def MaxShear(self, Direction):
        """
        Returns the maximum shear in the member for the given direction
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the maximum shear. Must be one of the following:
                "Fy" = Shear acting on the local y-axis
                "Fz" = Shear acting on the local z-axis
        """        
        
        Vmax = 0
        
        if Direction == "Fy":
            
            for segment in self.SegmentsZ:
                
                if segment.MaxShear() > Vmax:
                    
                    Vmax = segment.MaxShear()
                    
        if Direction == "Fz":
            
            for segment in self.SegmentsY:
                
                if segment.MaxShear() > Vmax:
                    
                    Vmax = segment.MaxShear()
        
        return Vmax
    
#%%
    def MinShear(self, Direction):
        """
        Returns the minimum shear in the member for the given direction
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the minimum shear. Must be one of the following:
                "Fy" = Shear acting on the local y-axis
                "Fz" = Shear acting on the local z-axis
        """        
        
        Vmin = 0
        
        if Direction == "Fy":
            
            for segment in self.SegmentsZ:
                
                if segment.MinShear() < Vmin:
                    
                    Vmin = segment.MinShear()
                    
        if Direction == "Fz":
            
            for segment in self.SegmentsY:
                
                if segment.MinShear() < Vmin:
                    
                    Vmin = segment.MinShear()
        
        return Vmin
    
#%%
    def PlotShear(self, Direction):
        """
        Plots the shear diagram for the member
        """
        
        ax = Member3D.__axis()
        ax.axhline(0, color='black', lw=1)
        ax.grid()
        
        x = []
        V = []
        
        # Calculate the shear diagram
        for i in range(20):
            x.append(self.L / 19 * i)
            V.append(self.Shear(Direction, self.L / 19 * i))

        Member3D.__plt.plot(x, V)
        Member3D.__plt.ylabel('Shear')
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.show()    
        
#%%
    def Moment(self, Direction, x):
        """
        Returns the moment at a point along the member's length
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the moment. Must be one of the following:
                "My" = Moment about the local y-axis
                "Mz" = moment about the local z-axis
        x : number
            The location at which to find the moment
        """
        
        # Check which axis is of interest
        if Direction == "My":
            
            # Check which segment "x" falls on
            for segment in self.SegmentsY:
                
                if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                    
                    return segment.Moment(x - segment.x1)
                
            if round(x,10) == round(self.L,10):
                
                lastIndex = len(self.SegmentsY) - 1
                return self.SegmentsY[lastIndex].Moment(x - self.SegmentsY[lastIndex].x1)
                
        elif Direction == "Mz":
            
            for segment in self.SegmentsZ:
                
                if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                    
                    return segment.Moment(x - segment.x1)
                
            if round(x,10) == round(self.L,10):
                
                lastIndex = len(self.SegmentsZ) - 1
                return self.SegmentsZ[lastIndex].Moment(x - self.SegmentsZ[lastIndex].x1)
            
#%%
    def MaxMoment(self, Direction):
        """
        Returns the maximum moment in the member for the given direction
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the maximum moment. Must be one of the following:
                "My" = Moment about the local y-axis
                "Mz" = Moment about the local z-axis
        """        
        
        Mmax = 0
        
        if Direction == "Mz":
            
            for segment in self.SegmentsZ:
                
                if segment.MaxMoment() > Mmax:
                    
                    Mmax = segment.MaxMoment()
                    
        if Direction == "My":
            
            for segment in self.SegmentsY:
                
                if segment.MaxMoment() > Mmax:
                    
                    Mmax = segment.MaxMoment()
        
        return Mmax
#%%
    def MinMoment(self, Direction):
        """
        Returns the minimum moment in the member for the given direction
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the minimum moment. Must be one of the following:
                "My" = Moment about the local y-axis
                "Mz" = Moment about the local z-axis
        """        
        
        Mmin = 0
        
        if Direction == "Mz":
            
            for segment in self.SegmentsZ:
                
                if segment.MinMoment() < Mmin:
                    
                    Mmin = segment.MinMoment()
                    
        if Direction == "My":
            
            for segment in self.SegmentsY:
                
                if segment.MinMoment() < Mmin:
                    
                    Mmin = segment.MinMoment()
        
        return Mmin
#%%
    def PlotMoment(self, Direction):
        """
        Plots the moment diagram for the member
        """
        
        ax = Member3D.__axis()
        ax.axhline(0, color='black', lw=1)
        ax.grid()
        
        x = []
        M = []
        
        # Calculate the moment diagram
        for i in range(20):
            
            x.append(self.L / 19 * i)
            M.append(self.Moment(Direction, self.L / 19 * i))

        Member3D.__plt.plot(x, M)
        Member3D.__plt.ylabel('Moment')
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.show()

#%%
    def Deflection(self, Direction, x):
        
        """
        Returns the deflection at a point along the member's length
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the deflection. Must be one of the following:
                "dy" = Deflection in the local y-axis
                "dz" = Deflection in the local x-axis
        x : number
            The location at which to find the deflection
        """
        
        # Check which axis is of interest
        if Direction == "dy":
            
            # Check which segment "x" falls on
            for segment in self.SegmentsZ:
                
                if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                    
                    return segment.Deflection(x - segment.x1)
                
            if round(x,10) == round(self.L,10):
                
                lastIndex = len(self.SegmentsZ) - 1
                return self.SegmentsZ[lastIndex].Deflection(x - self.SegmentsZ[lastIndex].x1)
                
        elif Direction == "dz":
            
            for segment in self.SegmentsY:
                
                if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                    
                    return segment.Deflection(x - segment.x1)
                
            if round(x,10) == round(self.L,10):
                
                lastIndex = len(self.SegmentsY) - 1
                return self.SegmentsY[lastIndex].Deflection(x - self.SegmentsY[lastIndex].x1) 

#%%
    def MaxDeflection(self, Direction):
        """
        Returns the maximum deflection in the member.
        
        Parameters
        ----------
        Direction : {"dy", "dz"}
            The direction in which to find the maximum deflection.
        """
        
        # Initialize the maximum deflection
        dmax = 0
        
        # Check the deflection at 100 locations along the member and find the largest value
        for i in range(100):
            d = self.Deflection(Direction, self.L * i / 99)
            if d > dmax:
                dmax = d
        
        # Return the largest value
        return dmax
    
#%%
    def MinDeflection(self, Direction):
        """
        Returns the minimum deflection in the member.
        
        Parameters
        ----------
        Direction : {"dy", "dz"}
            The direction in which to find the minimum deflection.
        """
        
        # Initialize the minimum deflection
        dmin = 0
        
        # Check the deflection at 100 locations along the member and find the smallest value
        for i in range(100):
            d = self.Deflection(Direction, self.L * i / 99)
            if d < dmin:
                dmin = d
        
        # Return the smallest value
        return dmin
              
#%%
    def PlotDeflection(self, Direction):
        """
        Plots the deflection diagram for the member
        """
        
        ax = Member3D.__axis()
        ax.axhline(0, color='black', lw=1)
        ax.grid()
        
        x = []
        d = []
        
        # Calculate the deflection diagram
        for i in range(20):
            
            x.append(self.L / 19 * i)
            d.append(self.Deflection(Direction, self.L / 19 * i))

        Member3D.__plt.plot(x, d)
        Member3D.__plt.ylabel('Deflection')
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.show()
        
#%%    
    # Divides the element up into mathematically continuous segments along each axis
    def SegmentMember(self):
        
        # Get the member's length and stiffness properties
        L = self.L
        E = self.M.E
        Iz = self.S.Izz
        Iy = self.S.Iyy
        SegmentsZ = self.SegmentsZ
        SegmentsY = self.SegmentsY
        
        # Create a list of discontinuity locations
        disconts = [0, L] # Member ends
        
        for load in self.PtLoads: 
            disconts.append(load[2]) # Point load locations
        
        for load in self.DistLoads: 
            disconts.append(load[3]) # Distributed load start locations
            disconts.append(load[4]) # Distributed load end locations
        
        # Sort the list and eliminate duplicate values
        disconts = sorted(set(disconts))
        
        # Clear out old data from any previous analyses
        SegmentsZ.clear()
        SegmentsY.clear()
        
        # Create a list of mathematically continuous segments for each direction
        for index in range(len(disconts) - 1):
            
            # z-direction segments (bending about local z-axis)
            newSeg = BeamSegment()        # Create the new segment
            newSeg.x1 = disconts[index]   # Segment start location
            newSeg.x2 = disconts[index+1] # Segment end location
            newSeg.EI = E*Iz              # Segment flexural stiffness
            SegmentsZ.append(newSeg)      # Add the segment to the list
            
            # y-direction segments (bending about local y-axis)
            newSeg = BeamSegment()        # Create the new segment
            newSeg.x1 = disconts[index]   # Segment start location
            newSeg.x2 = disconts[index+1] # Segment end location
            newSeg.EI = E*Iy              # Segment flexural stiffness
            SegmentsY.append(newSeg)      # Add the segment to the list
        
        # Get the member local end forces, local fixed end reactions, and local displacements
        f = self.f()     # Member local end force vector
        fer = self.__fer_Unc() # Member local fixed end reaction vector
        d = self.d()     # Member local displacement vector
        
        # Get the local deflections and calculate the slope at the start of the member
        # Note that the slope may not be available directly from the local displacement vector if member end releases have been used,
        # so slope-deflection has been applied to solve for it
        m1z = f[5, 0]       # local z-axis moment at start of member
        m2z = f[11, 0]      # local z-axis moment at end of member
        m1y = -f[4, 0]       # local y-axis moment at start of member
        m2y = -f[10, 0]      # local y-axis moment at end of member
        fem1z = fer[5, 0]   # local z-axis fixed end moment at start of member
        fem2z = fer[11, 0]  # local z-axis fixed end moment at end of member
        fem1y = -fer[4, 0]   # local y-axis fixed end moment at start of member
        fem2y = -fer[10, 0]  # local y-axis fixed end moment at end of member
        delta1y = d[1, 0]   # local y displacement at start of member
        delta2y = d[7, 0]   # local y displacement at end of member
        delta1z = d[2, 0]   # local z displacement at start of member
        delta2z = d[8, 0]   # local z displacement at end of member
        SegmentsZ[0].delta1 = delta1y
        SegmentsY[0].delta1 = delta1z
        SegmentsZ[0].theta1 = 1/3*((m1z - fem1z)*L/(E*Iz) - (m2z - fem2z)*L/(2*E*Iz) + 3*(delta2y - delta1y)/L)
        SegmentsY[0].theta1 = 1/3*((m1y - fem1y)*L/(E*Iy) - (m2y - fem2y)*L/(2*E*Iy) + 3*(delta2z - delta1z)/L)
        
        # Add loads to each segment
        for i in range(len(SegmentsZ)):
            
            # Get the starting point of the segment
            x = SegmentsZ[i].x1
            
            # Initialize the distributed loads on the segment to zero
            SegmentsZ[i].w1 = 0
            SegmentsZ[i].w2 = 0
            SegmentsZ[i].p1 = 0
            SegmentsZ[i].p2 = 0
            SegmentsY[i].w1 = 0
            SegmentsY[i].w2 = 0
            SegmentsY[i].p1 = 0
            SegmentsY[i].p2 = 0
            
            # Initialize the slope and displacement at the start of the segment
            if i > 0: # The first segment has already been initialized
                SegmentsZ[i].theta1 = SegmentsZ[i-1].Slope(SegmentsZ[i-1].Length())
                SegmentsZ[i].delta1 = SegmentsZ[i-1].Deflection(SegmentsZ[i-1].Length())
                SegmentsY[i].theta1 = SegmentsY[i-1].Slope(SegmentsY[i-1].Length())
                SegmentsY[i].delta1 = SegmentsY[i-1].Deflection(SegmentsY[i-1].Length())
                
            # Add the effects of the beam end forces to the segment
            SegmentsZ[i].P1 = f[0,0]
            SegmentsZ[i].V1 = f[1,0]
            SegmentsZ[i].M1 = -f[5,0] + f[1,0]*x
            SegmentsY[i].P1 = f[0,0]
            SegmentsY[i].V1 = f[2,0]
            SegmentsY[i].M1 = f[4,0] + f[2,0]*x
            
            # Add effects of point loads occuring prior to this segment
            for ptLoad in self.PtLoads:
                
                if round(ptLoad[2],10) <= round(x,10):
                    
                    if ptLoad[0] == "Fx":
                        SegmentsZ[i].P1 += ptLoad[1]
                    elif ptLoad[0] == "Fy":
                        SegmentsZ[i].V1 -= ptLoad[1]
                        SegmentsZ[i].M1 -= ptLoad[1]*(x-ptLoad[2])
                    elif ptLoad[0] == "Fz":
                        SegmentsY[i].V1 -= ptLoad[1]
                        SegmentsY[i].M1 -= ptLoad[1]*(x - ptLoad[2])
                    elif ptLoad[0] == "My":
                        SegmentsY[i].M1 -= ptLoad[1]
                    elif ptLoad[0] == "Mz":
                        SegmentsZ[i].M1 -= ptLoad[1]
            
            # Add distributed loads to the segment
            for distLoad in self.DistLoads:
                
                # Get the parameters for the distributed load
                Direction = distLoad[0]
                w1 = distLoad[1]
                w2 = distLoad[2]
                x1 = distLoad[3]
                x2 = distLoad[4]
            
                # Determine if the load affects the segment
                if round(x1,10) <= round(x,10):
                    
                    if Direction == "Fx":
                        
                        # Determine if the load ends after the start of the segment
                        if round(x2,10) > round(x,10):
                                                
                            # Break up the load and place it on the segment
                            # Note that 'w1' and 'w2' are really 'p1' and 'p2' here
                            SegmentsZ[i].p1 += (w2 - w1)/(x2 - x1)*(x - x1) + w1
                            SegmentsZ[i].p2 += (w2 - w1)/(x2 - x1)*(SegmentsZ[i].x2 - x1) + w1
                            SegmentsY[i].p1 += (w2 - w1)/(x2 - x1)*(x - x1) + w1
                            SegmentsY[i].p2 += (w2 - w1)/(x2 - x1)*(SegmentsY[i].x2 - x1) + w1
                            
                            # Calculate the magnitude of the load at the start of the segment
                            w2 = w1+(w2-w1)/(x2-x1)*(x-x1)
                            x2 = x
                        
                        # Calculate the axial force at the start of the segment
                        SegmentsZ[i].P1 -= (w1 + w2)/2*(x2 - x1)
                        SegmentsY[i].P1 -= (w1 + w2)/2*(x2 - x1)
                    
                    elif Direction == "Fy":
                        
                        # Determine if the load ends after the start of the segment
                        if round(x2,10) > round(x,10):
                                                
                            # Break up the load and place it on the segment
                            SegmentsZ[i].w1 += (w2 - w1)/(x2 - x1)*(x - x1) + w1
                            SegmentsZ[i].w2 += (w2 - w1)/(x2 - x1)*(SegmentsZ[i].x2 - x1) + w1
                            
                            # Calculate the magnitude of the load at the start of the segment
                            # This will be used as the 'x2' value for the load prior to the start of the segment
                            w2 = w1 + (w2 - w1)/(x2 - x1)*(x - x1)
                            x2 = x
                        
                        # Calculate the shear and moment at the start of the segment due to the load
                        SegmentsZ[i].V1 -= (w1 + w2)/2*(x2 - x1)
                        SegmentsZ[i].M1 -= (x1 - x2)*(2*w1*x1 - 3*w1*x + w1*x2 + w2*x1 - 3*w2*x + 2*w2*x2)/6
                    
                    elif Direction == "Fz":
                        
                        # Determine if the load ends after the start of the segment
                        if round(x2,10) > round(x,10):
                                                
                            # Break up the load and place it on the segment
                            SegmentsY[i].w1 += (w2 - w1)/(x2 - x1)*(SegmentsY[i].x1 - x1) + w1
                            SegmentsY[i].w2 += (w2 - w1)/(x2 - x1)*(SegmentsY[i].x2 - x1) + w1
                            
                            # Calculate the magnitude of the load at the start of the segment
                            w2 = w1 + (w2 - w1)/(x2 - x1)*(x - x1)
                            x2 = x
                        
                        # Calculate the shear and moment at the start of the segment due to the load
                        SegmentsY[i].V1 -= (w1 + w2)/2*(x2 - x1)
                        SegmentsY[i].M1 -= (x1 - x2)*(2*w1*x1 - 3*w1*x + w1*x2 + w2*x1 - 3*w2*x + 2*w2*x2)/6
