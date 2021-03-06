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

parser.add_argument('--draw-stress',  help='Draw 3D frame of members indicating stress.', default=None,          choices=['seq', 'sxx', 'txy', 'tzx'])
parser.add_argument('--draw-frame',   help='Draw 3D frame of members.',                                          action='store_true')
parser.add_argument('--wire-frame',   help='Draw 3D wire frame of members.',                                     action='store_true')
parser.add_argument('--plot-results', help='Plot shear, moment & displacement results (left; vertical).',        action='store_true')
parser.add_argument('--lhs-rotate',   help='Rotate beam through 90 degrees (left).',                             action='store_true')
parser.add_argument('--rhs-rotate',   help='Rotate beam through 90 degrees (right).',                            action='store_true')
parser.add_argument('--lhs-section',  help='Specify section for beam on left  [RHS].',    default='RHS',         choices=['Rectangle', 'RHS', 'Circular', 'CHS', 'Universal', 'BB.RHS.Cold', 'BB.CHS.Hot', 'BB.UB'])
parser.add_argument('--rhs-section',  help='Specify section for beam on right [RHS].',    default='RHS',         choices=['Rectangle', 'RHS', 'Circular', 'CHS', 'Universal', 'BB.RHS.Cold', 'BB.CHS.Hot', 'BB.UB'])
parser.add_argument('--lhs-rounded',  help='In case of Rectangular/RHS section, round the corners.',             action='store_true')
parser.add_argument('--rhs-rounded',  help='In case of Rectangular/RHS section, round the corners.',             action='store_true')
parser.add_argument('--lhs-material', help='Specify material for beam on left  [steel].', default='steel',       choices=['steel', 'aluminium'])
parser.add_argument('--rhs-material', help='Specify material for beam on right [steel].', default='steel',       choices=['steel', 'aluminium'])
parser.add_argument('--case',         help='Specify load case [U2B.Z.-100k].',            default='U2B.Z.-100k', choices=['U2B.Z.-100k','U2B.Y.-100k','U4B.Z.-+50k','U4B.Y.-+50k','F4B.Z.+-10k','F4B.Y.+-10k'])

args = parser.parse_args()

if args.case == 'U2B.Z.-100k':
    case_TwoBeams = True
    case_Y_Load   = 0
    case_Z_Load   = -100E3
    case_EndFixed = False
elif args.case == 'U2B.Y.-100k':
    case_TwoBeams = True
    case_Y_Load   = -100E3
    case_Z_Load   = 0
    case_EndFixed = False
elif args.case == 'U4B.Z.-+50k':
    case_TwoBeams = False
    case_Y_Load   = 0
    case_Z_Load   = -50E3
    case_EndFixed = False
elif args.case == 'U4B.Y.-+50k':
    case_TwoBeams = False
    case_Y_Load   = -50E3
    case_Z_Load   = 0
    case_EndFixed = False
elif args.case == 'F4B.Z.+-10k':
    case_TwoBeams = False
    case_Y_Load   = 0
    case_Z_Load   = 10E3
    case_EndFixed = True
else: # args.case == 'F4B.Y.+-10k':
    case_TwoBeams = False
    case_Y_Load   = 10E3
    case_Z_Load   = 0
    case_EndFixed = True

# Option to twist beam through 90 degrees
# default: unless the beam itself is vertical, the beam and gravity define the zx-plane
if args.lhs_rotate:
    lhs_ref = [0,10,0.1]
else:
    lhs_ref = None

if args.rhs_rotate:
    rhs_ref = [0,10,0.1]
else:
    rhs_ref = None

UB    = Section.BlueBook.UniversalBeam()
#CHS_C = Section.BlueBook.CHS_ColdFormed()
CHS_H = Section.BlueBook.CHS_HotFinished()
RHS_C = Section.BlueBook.RHS_ColdFormed()
#RHS_H = Section.BlueBook.RHS_HotFinished()
#SHS_C = Section.BlueBook.SHS_ColdFormed()
#SHS_H = Section.BlueBook.SHS_HotFinished()

if args.lhs_section == 'Rectangle':
    if args.lhs_rounded:
        lhs_section = Section.Rectangular(0.3, 0.1, 0.01) # 300mm x 100mm x 10mm
    else:
        lhs_section = Section.Rectangular(0.3, 0.1) # 300mm x 100mm
elif args.lhs_section == 'RHS':
    if args.lhs_rounded:
        lhs_section = Section.RHS(0.1, 0.1, 0.005, 0.01) # 100mm x 100mm x 10mm, 5mm thick
    else:
        lhs_section = Section.RHS(0.1, 0.1, 0.005) # 100mm x 100mm, 5mm thick
elif args.lhs_section == 'Circular':
    lhs_section = Section.Circular(0.03) # 30mm
elif args.lhs_section == 'CHS':
    lhs_section = Section.CHS(0.1, 0.005) # 100mm diameter, 5mm thick
elif args.lhs_section == 'BB.RHS.Cold':
    lhs_section = RHS_C['100  x  50']['5.0']
elif args.lhs_section == 'BB.CHS.Hot':
    lhs_section = CHS_H['76.1']['5.0']
elif args.lhs_section == 'BB.UB':
    lhs_section = UB['203 x 133']['x 30']
else: # 'Universal'
    lhs_section = Section.Universal(0.1138, 0.2068, 0.0096, 0.0063, 0.0076) # breadth / depth / flange / web [/ root radius]

if args.rhs_section == 'Rectangle':
    if args.rhs_rounded:
        rhs_section = Section.Rectangular(0.3, 0.1, 0.01) # 300mm x 100mm x 10mm
    else:
        rhs_section = Section.Rectangular(0.3, 0.1) # 300mm x 100mm
elif args.rhs_section == 'RHS':
    if args.rhs_rounded:
        rhs_section = Section.RHS(0.1, 0.1, 0.005, 0.01) # 100mm x 100mm x 10mm, 5mm thick
    else:
        rhs_section = Section.RHS(0.1, 0.1, 0.005) # 100mm x 100mm, 5mm thick
elif args.rhs_section == 'Circular':
    rhs_section = Section.Circular(0.03) # 30mm
elif args.rhs_section == 'CHS':
    rhs_section = Section.CHS(0.1, 0.005) # 100mm diameter, 5mm thick
elif args.rhs_section == 'BB.RHS.Cold':
    rhs_section = RHS_C['100  x  50']['5.0']
elif args.rhs_section == 'BB.CHS.Hot':
    rhs_section = CHS_H['76.1']['5.0']
elif args.rhs_section == 'BB.UB':
    rhs_section = UB['203 x 133']['x 30']
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

# Add nodes (for total beam length of 3 metres along x-axis, but offset to avoid origin)
L = 3
SimpleBeam.AddNode("N1", -L/2, 0.1, 0.1)
SimpleBeam.AddNode("N2",  0,   0.1, 0.1)
SimpleBeam.AddNode("N3",  L/2, 0.1, 0.1)

if not case_TwoBeams:
    SimpleBeam.AddNode("N4", -L/4, 0.1, 0.1)
    SimpleBeam.AddNode("N5",  L/4, 0.1, 0.1)

    SimpleBeam.AddMemberExt("M1", "N1", "N4", lhs_material, lhs_section, lhs_ref) # *_ref: optional argument
    SimpleBeam.AddMemberExt("M2", "N4", "N2", lhs_material, lhs_section, lhs_ref) # *_ref: optional argument
    SimpleBeam.AddMemberExt("M3", "N2", "N5", rhs_material, rhs_section, rhs_ref) #        [default: None]
    SimpleBeam.AddMemberExt("M4", "N5", "N3", rhs_material, rhs_section, rhs_ref) #        [default: None]

    # Add opposite point loads at the 1/4 & 3/4 points of the beam
    if case_Y_Load != 0:
        SimpleBeam.AddNodeLoad("N4", "FY",  case_Y_Load)
        SimpleBeam.AddNodeLoad("N5", "FY", -case_Y_Load)
    if case_Z_Load != 0:
        SimpleBeam.AddNodeLoad("N4", "FZ",  case_Z_Load)
        SimpleBeam.AddNodeLoad("N5", "FZ", -case_Z_Load)
else:
    SimpleBeam.AddMemberExt("M1", "N1", "N2", lhs_material, lhs_section, lhs_ref) # *_ref: optional argument
    SimpleBeam.AddMemberExt("M2", "N2", "N3", rhs_material, rhs_section, rhs_ref) #        [default: None]

    # Add a point load at the midspan of the beam
    if case_Y_Load != 0:
        SimpleBeam.AddNodeLoad("N2", "FY", case_Y_Load)
    if case_Z_Load != 0:
        SimpleBeam.AddNodeLoad("N2", "FZ", case_Z_Load)

if case_EndFixed:
    # Provide simple supports - completely fixed ends
    SimpleBeam.DefineSupport("N1", True, True, True, True, True, True)
    SimpleBeam.DefineSupport("N3", True, True, True, True, True, True)
else:
    # Provide simple supports - can't displace or twist at either end
    SimpleBeam.DefineSupport("N1", True, True, True, True, False, False)
    SimpleBeam.DefineSupport("N3", True, True, True, True, False, False)

if args.wire_frame:
    # Draw interactive wire-frame showing sections
    SimpleBeam.Display(True)
elif args.draw_frame:
    # Draw interactive frame showing sections
    SimpleBeam.Display(False)

if case_TwoBeams:
    # Expected moment
    MY = -case_Z_Load * L / 4
    MZ = -case_Y_Load * L / 4
    # Expected deflection - only if materials and sections are same and same orientation
    if args.lhs_rotate:
        delta_Y = case_Y_Load * L**3 / (48 * lhs_material.E * lhs_section.Iyy)
        delta_Z = case_Z_Load * L**3 / (48 * lhs_material.E * lhs_section.Izz)
    else:
        delta_Y = case_Y_Load * L**3 / (48 * lhs_material.E * lhs_section.Izz)
        delta_Z = case_Z_Load * L**3 / (48 * lhs_material.E * lhs_section.Iyy)
    print('Expected deflection: dy = {dy:.2f} mm, dz = {dz:.2f} mm; max. bending moment: My = {my:.2f} kNm, Mz = {mz:.2f} kNm.'.format(dy=(delta_Y * 1E3), dz=(delta_Z * 1E3), my=(MY / 1E3), mz=(MZ / 1E3)))

# Analyze the beam
SimpleBeam.Analyze()
print('Analysed deflection: dy = {dy:.2f} mm, dz = {dz:.2f} mm.'.format(dy=(SimpleBeam.GetNode("N2").DY * 1E3), dz=(SimpleBeam.GetNode("N2").DZ * 1E3)))

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
print("Left  Support Reaction: ({ry:.2f}, {rz:.2f}) kN".format(ry=(SimpleBeam.GetNode("N1").RxnFY / 1E3), rz=(SimpleBeam.GetNode("N1").RxnFZ / 1E3)))
print("Right Support Reaction: ({ry:.2f}, {rz:.2f}) kN".format(ry=(SimpleBeam.GetNode("N3").RxnFY / 1E3), rz=(SimpleBeam.GetNode("N3").RxnFZ / 1E3)))

if args.draw_stress is not None:
    # Draw interactive wire-frame showing sections
    SimpleBeam.DisplayResults(args.draw_stress)
