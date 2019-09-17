# -*- coding: utf-8 -*-
"""
Created on Sun Sep 16 14:43:31 2018
Modified on Mon Sep 17 2019
@author: craig, fjf
"""
# Example of a simply supported beam with a point load.
# Units used in this example are metric (SI)
# 
# To run, e.g.:
# 
# PYTHONPATH=$PWD python3 Examples/Double\ Beam\ -\ Point\ Load.py --help
# 
# Requirements (via pip modules):
#   PyNite         NumPy
#   --draw-frame   VisPy
#   --plot-results MatPlotLib


import argparse

import numpy as np

from PyNite import FEModel3D, Section
from PyNite.Material import Material

parser = argparse.ArgumentParser(description="Beam simply supported at either end, loaded in the middle.")

parser.add_argument('--draw-frame',   help='Draw 3D wire frame of members.',                               action='store_true')
parser.add_argument('--plot-results', help='Plot shear, moment & displacement results (left; vertical).',  action='store_true')
parser.add_argument('--lhs-rotate',   help='Rotate beam through 90 degrees (left).',                       action='store_true')
parser.add_argument('--rhs-rotate',   help='Rotate beam through 90 degrees (right).',                      action='store_true')
parser.add_argument('--lhs-section',  help='Specify section for beam on left  [RHS].',    default='RHS',   choices=['Rectangle', 'Rounded', 'RHS', 'Circular', 'CHS', 'Universal'])
parser.add_argument('--rhs-section',  help='Specify section for beam on right [RHS].',    default='RHS',   choices=['Rectangle', 'Rounded', 'RHS', 'Circular', 'CHS', 'Universal'])
parser.add_argument('--lhs-material', help='Specify material for beam on left  [steel].', default='steel', choices=['steel', 'aluminium'])
parser.add_argument('--rhs-material', help='Specify material for beam on right [steel].', default='steel', choices=['steel', 'aluminium'])

args = parser.parse_args()

# Option to twist beam through 90 degrees
# default: unless the beam itself is vertical, the beam and gravity define the zx-plane
if args.lhs_rotate:
    lhs_ref = [0,10,0]
else:
    lhs_ref = None

if args.rhs_rotate:
    rhs_ref = [0,10,0]
else:
    rhs_ref = None

if args.lhs_section == 'Rectangle':
    lhs_section = Section.Rectangular(0.3, 0.1) # 300mm x 100mm
elif args.lhs_section == 'Rounded':
    lhs_section = Section.Rectangular(0.2, 0.2, 0.01) # 200mm x 200mm x 10mm
elif args.lhs_section == 'RHS':
    lhs_section = Section.RHS(0.3, 0.1, 0.005, 0.01) # 300mm x 100mm x 10mm, 5mm thick
elif args.lhs_section == 'Circular':
    lhs_section = Section.Circular(0.4) # 400mm
elif args.lhs_section == 'CHS':
    lhs_section = Section.CHS(0.5, 0.03) # 500mm, 30mm thick
else: # 'Universal'
    lhs_section = Section.Universal(0.1138, 0.2068, 0.0096, 0.0063, 0.0076) # breadth / depth / flange / web [/ root radius]

if args.rhs_section == 'Rectangle':
    rhs_section = Section.Rectangular(0.3, 0.1) # 300mm x 100mm
elif args.rhs_section == 'Rounded':
    rhs_section = Section.Rectangular(0.2, 0.2, 0.01) # 200mm x 200mm x 10mm
elif args.rhs_section == 'RHS':
    rhs_section = Section.RHS(0.3, 0.1, 0.005, 0.01) # 300mm x 100mm x 10mm, 5mm thick
elif args.rhs_section == 'Circular':
    rhs_section = Section.Circular(0.4) # 400mm
elif args.rhs_section == 'CHS':
    rhs_section = Section.CHS(0.5, 0.03) # 500mm, 30mm thick
else: # 'Universal'
    rhs_section = Section.Universal(0.1138, 0.2068, 0.0096, 0.0063, 0.0076) # breadth / depth / flange / web [/ root radius]

if args.lhs_material == 'steel':
    lhs_material = Material.steel()
else: # 'aluminium'
    lhs_material = Material.aluminium()

if args.rhs_material == 'steel':
    rhs_material = Material.steel()
else: # 'aluminium'
    rhs_material = Material.aluminium()

# Create a new finite element model
SimpleBeam = FEModel3D()

# Add nodes (3 metres apart on x-axis)
L = 3
SimpleBeam.AddNode("N1", -L/2, 0.1, 0.1)
SimpleBeam.AddNode("N2",  0,   0.1, 0.1)
SimpleBeam.AddNode("N3",  L/2, 0.1, 0.1)

SimpleBeam.AddMemberExt("M1", "N1", "N2", lhs_material, lhs_section, lhs_ref) # *_ref: optional argument
SimpleBeam.AddMemberExt("M2", "N2", "N3", rhs_material, rhs_section, rhs_ref) #        [default: None]

# Provide simple supports - can't displace or twist at either end
SimpleBeam.DefineSupport("N1", True, True, True, True, False, False)
SimpleBeam.DefineSupport("N3", True, True, True, True, False, False)

# Add a point load at the midspan of the beam
F = -1E3 # load, e.g., a hanging weight of approx 100 kg
SimpleBeam.AddNodeLoad("N2", "FZ", F)

if args.draw_frame:
    # Draw interactive wire-frame showing sections
    SimpleBeam.Display()

# Expected moment
M = -F * L / 4
# Expected deflection - only if materials and sections are same and same orientation
if args.lhs_rotate:
    delta = F * L**3 / (48 * lhs_material.E * lhs_section.Izz)
else:
    delta = F * L**3 / (48 * lhs_material.E * lhs_section.Iyy)
print('Expected deflection = {d:.2f} mm; max. bending moment = {m:.2f}'.format(d=(delta * 1E3), m=M))

# Analyze the beam
SimpleBeam.Analyze()
print('Deflection from analysis = {d:.2f} mm'.format(d=(SimpleBeam.GetNode("N2").DZ * 1E3)))

if args.plot_results:
    # Plot the shear, moment, and deflection diagrams
    if args.lhs_rotate:
        SimpleBeam.GetMember("M1").PlotShear("Fy")
        SimpleBeam.GetMember("M1").PlotMoment("Mz")
        SimpleBeam.GetMember("M1").PlotDeflection("dy")
    else:
        SimpleBeam.GetMember("M1").PlotShear("Fz")
        SimpleBeam.GetMember("M1").PlotMoment("My")
        SimpleBeam.GetMember("M1").PlotDeflection("dz")

# Print reactions at each end of the beam
print("Left Support Reaction: {Rxn:.2f} N".format(Rxn = SimpleBeam.GetNode("N1").RxnFZ))
print("Right Support Reacton: {Rxn:.2f} N".format(Rxn = SimpleBeam.GetNode("N3").RxnFZ))
