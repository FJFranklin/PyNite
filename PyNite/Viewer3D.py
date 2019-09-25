# -*- coding: utf-8 -*-
"""
Created September 2019
@author: fjf
"""
# %%
import numpy as np

# %%
class Viewer3D(object):
    """
    A class for displaying a 3D model.
    """
#%%
    scaling = 10    # VisPy rendering very poor when values < 1, so scale up
    __scene = None  # Don't import VisPy unless and until actually needed
    __geometry = None
#%% 
    def __init__(self):
        if Viewer3D.__scene is None:
            from vispy import scene, geometry
            Viewer3D.__scene     = scene
            Viewer3D.__geometry  = geometry

            # Available colormaps in Vispy:
            # autumn, blues, cool, greens, reds, spring, summer, fire, grays, hot, ice, winter, light_blues,
            # orange, viridis, coolwarm, PuGr, GrBu, GrBu_d, RdBu, cubehelix, single_hue, hsl, husl, diverging, RdYeBuCy

            from vispy.color import Colormap
            self.__cmap = Colormap(['blue','cyan','green','yellow','red'])

        self.__canvas = Viewer3D.__scene.SceneCanvas(keys='interactive', size=(1200, 900), show=True)

        # Set up a viewbox to display the cube with interactive arcball
        self.__view = self.__canvas.central_widget.add_view()
        self.__view.bgcolor = '#cfcfef'
        self.__view.camera  = Viewer3D.__scene.cameras.turntable.TurntableCamera(scale_factor=10)
        self.__view.padding = 10

        Viewer3D.__scene.visuals.XYZAxis(parent=self.__view.scene)

        self.__vertex = None
        self.__vcolor = None

#%% 
    def Run(self, smooth=True):
        if self.__vertex is not None:
            data = Viewer3D.__geometry.MeshData(vertices=self.__vertex, vertex_colors=self.__vcolor)
            if smooth:
                mesh = Viewer3D.__scene.visuals.Mesh(meshdata=data, shading='smooth', parent=self.__view.scene)
            else:
                mesh = Viewer3D.__scene.visuals.Mesh(meshdata=data, parent=self.__view.scene)

        self.__canvas.app.run()

#%% 
    def Line(self, vertices, color=None):
        """
        Draws a line in 3D

        Parameters
        ----------
        vertices : numpy.ndarray
            Array of vertices (Nv,3) to draw line through
        color : [r,g,b,a] array
            Color to use for drawing line; optional [default: black]
        """
        if color is None:
            color = 'black'

        Viewer3D.__scene.visuals.Line(pos=(vertices * Viewer3D.scaling), color=color, parent=self.__view.scene)

#%% 
    def Add(self, triangle, color, uniform=True):
        """
        Adds a triangle to the 3D mesh

        Parameters
        ----------
        triangle : numpy.ndarray
            Array of vertices (3,3) for triangle
        color : [[r,g,b,a],] array
            Array of one or three colors to use for displaying face
        """
        if self.__vertex is None:
            self.__vertex = np.asarray([triangle]) * Viewer3D.scaling
            if len(color) == 1:
                self.__vcolor = np.zeros((1,3,4))
                self.__vcolor[0,:] = color[0]
            else:
                self.__vcolor = np.asarray([color])
        else:
            vcoord = np.asarray([triangle]) * Viewer3D.scaling
            if len(color) == 1:
                vcolor = np.zeros((1,3,4))
                vcolor[0,:] = color[0]
            else:
                vcolor = np.asarray([color])
            self.__vertex = np.append(self.__vertex, vcoord, axis=0)
            self.__vcolor = np.append(self.__vcolor, vcolor, axis=0)

#%% 
    def ColorFromValue(self, c):
        """
        Converts a value in a range to a color

        Parameters
        ----------
        c : numpy.ndarray
            Array of values to map to color

        Returns
        -------
        color : [[r,g,b,a],] array
            mapped color
        """
        return self.__cmap.map(c)

#%% 
    def ColorBar(self, description, c_min_max_str):
        """
        Add a colorbar to the plot

        """
        # Overlays a grid to position subviews & widgets:
        grid = self.__canvas.central_widget.add_grid(margin=10)

        # Positioning of text in colorbar is unhelpful; this hackery seems to make it look okay, however
        cbar = Viewer3D.__scene.widgets.ColorBarWidget(self.__cmap, 'left', clim=c_min_max_str, label=('\n\n '+description), axis_ratio=0.075, padding=[0.5,0.5], border_width=1)

        grid.add_widget(cbar, col=0)
        # which centers the colorbar unless we add something else to effectively nudge it to the left:
        view = grid.add_view(col=1, col_span=3)
        view.visible = False
