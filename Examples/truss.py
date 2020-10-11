import numpy as np

# Materials: PyNite has steel and aluminium predefined, and you can make your own

from PyNite.Material import Material

def make_gold():
    gold = Material(77.2E9, 27.2E9, 19320) # Young's elastic modulus; shear modulus; density
    gold.set_color((0.9, 0.9, 0.0))        # yellow - for display purposes only
    gold.sigma_y = 47E6                    # estimated yield stress - very alloy-dependent
    return gold

material_1 = make_gold()
material_2 = Material.steel()     # Defaults provided for steel and aluminium
material_3 = Material.aluminium() # (Only isotropic elastic materials supported)

# Cross-sections: You can custom-define solid and hollow rectangular and circular sections, and I-beams

from PyNite import Section

def rectangular_solid_section(width, height, corner_radius=0): # [units: mm]
    return Section.Rectangular(width/1E3, height/1E3, corner_radius/1E3)

def rectangular_hollow_section(width, height, wall_thickness, corner_radius=0): # [units: mm]
    return Section.RHS(width/1E3, height/1E3, wall_thickness/1E3, corner_radius/1E3)

def circular_solid_section(diameter): # [units: mm]
    return Section.Circular(diameter/1E3)

def circular_hollow_section(diameter, wall_thickness): # [units: mm]
    return Section.CHS(diameter/1E3, wall_thickness/1E3)

def universal_section(width, height, flange_thickness, web_thickness, corner_radius=0): # I-beam [units: mm]
    return Section.Universal(width/1E3, height/1E3, flange_thickness/1E3, web_thickness/1E3, corner_radius/1E3)

section_1 = rectangular_solid_section(10, 40, 3)
section_2 = rectangular_hollow_section(40, 40, 3, 5)
section_3 = circular_solid_section(30)
section_4 = circular_hollow_section(50, 2.5)
section_5 = universal_section(40, 60, 3, 2, 1)

# ... or use one from sets of standard industrial sections defined in the steel industry's Blue Book

#CHS_C = Section.BlueBook.CHS_ColdFormed()                                                                                                                                                                                      
#CHS_H = Section.BlueBook.CHS_HotFinished()
#RHS_C = Section.BlueBook.RHS_ColdFormed()
#RHS_H = Section.BlueBook.RHS_HotFinished()                                                                                                                                                                                     
#SHS_C = Section.BlueBook.SHS_ColdFormed()                                                                                                                                                                                      
#SHS_H = Section.BlueBook.SHS_HotFinished()                                                                                                                                                                                     
UB    = Section.BlueBook.UniversalBeam()

section_prop_header = '|Iyy [cm4]|Izz [cm4]| J [cm4]|Area [cm2]'
section_prop_divide = '|---------|---------|--------|----------'

def section_properties(section):
    Iyy = section.Iyy * 1E8
    Izz = section.Izz * 1E8
    J   = section.J   * 1E8
    A   = section.A   * 1E4
    return '| {y:>7.0f} | {z:>7.1f} | {J:>6.1f} |  {A:>6.1f}  '.format(y=Iyy, z=Izz, J=J, A=A)

def section_props_list(section_list):
    print(' # ' + section_prop_header)
    print('---' + section_prop_divide)
    s_no = 0
    for s in section_list:
        s_no += 1
        print('{s:>2d} '.format(s=s_no) + section_properties(s))
    print('---' + section_prop_divide)
    print() # blank line

def UB_list():
    page = Section.BlueBook.UniversalBeam()
    print('Designation |Mass / m|Width [mm]|Height [mm]|Flange [mm]|Web [mm]|Radius [mm]' + section_prop_header)
    print('------------|--------|----------|-----------|-----------|--------|-----------' + section_prop_divide)
    for d in page.keys():
        designation = d
        sections = page[d]
        for s in sections.keys():
            section = sections[s]
            width  = section.breadth * 1E3 # convert to millimetres
            height = section.depth   * 1E3
            flange = section.flange  * 1E3
            web    = section.web     * 1E3
            radius = section.radius  * 1E3
            fmtstr = ' {d:>10s} | {m:6} |  {w:>6.1f}  |  {h:>7.1f}  |  {f:>7.1f}  |  {b:>4.1f}  |  {r:>7.1f}  '
            dimstr = fmtstr.format(d=designation, m=s, w=width, h=height, f=flange, b=web, r=radius)
            props = section_properties(section)
            print(dimstr + props)
            designation = ''
        print('------------|--------|----------|-----------|-----------|--------|-----------' + section_prop_divide)
    print() # blank line

#To list the available sections in UB:
#UB_list()

section_6 = UB['127 x 76']['x 13']

section_props_list([section_1, section_2, section_3, section_4, section_5, section_6])

from PyNite import FEModel3D

def accelerate_member(model, member_name, acceleration=None):
    if acceleration is None:
        acceleration = [0,0,9.81] # default is to add weight; gravity is equivalent to upward acceleration

    member = model.GetMember(member_name)

    material = member.M
    section  = member.S
    length   = member.L
    area     = section.A
    density  = material.density

    inertial = -np.asarray(acceleration) * density * area

    wx = np.dot(inertial, member.e_i) # transform applied force to beam-local coordinate system
    wy = np.dot(inertial, member.e_j)
    wz = np.dot(inertial, member.e_k)

    model.AddMemberDistLoad(member_name, 'Fx', -wx, -wx, 0, length) # uniform along the length; sign convention
    model.AddMemberDistLoad(member_name, 'Fy', -wy, -wy, 0, length)
    model.AddMemberDistLoad(member_name, 'Fz', -wz, -wz, 0, length)
    
    return np.asarray([wx, wy, wz])

def cantilever(section, material, length, force, moment, **kwargs):
    """Analyse a cantilever bending problem
    """
    # Add some options to display frame & results
    if 'draw' in kwargs: # Note: requires VisPy [options: 'draw', 'frame']
        draw = kwargs['draw']
    else:
        draw = None
    if 'plot' in kwargs: # [options: True, False]
        plot = kwargs['plot']
    else:
        plot = False
    if 'stress' in kwargs: # Note: experimental; requires VisPy [options: 'seq', 'sxx', 'txy', 'tzx']
        stress = kwargs['stress']
    else:
        stress = None
    if 'gravity' in kwargs: # [options: True, False]
        gravity = kwargs['gravity']
    else:
        gravity = False

    # Create a new finite element model                                                                                                                                                                                             
    model = FEModel3D()

    # Add nodes (N1 and N2) at ends of beam; offset from origin for the sake of visualisation                                                                                                                                        
    model.AddNode('N1', 0,      0.1, 0.1)
    model.AddNode('N2', length, 0.1, 0.1)

    # Add a beam (M1) between the nodes
    model.AddMemberExt('M1', 'N1', 'N2', material, section)

    # Fix to wall - can't move, can't twist, zero gradient
    model.DefineSupport('N1', True, True, True, True, True, True)

    if force != 0:
        # Add a vertical (Z) load at the end
        model.AddNodeLoad('N2', 'FZ', force)

    if moment != 0:
        # Add a moment at the end in the ZX-plane
        model.AddNodeLoad('N2', 'MY', moment)

    if gravity:
        w = accelerate_member(model, 'M1')
        uniform = w[2] # z-component of weight per unit length
    else:
        uniform = 0

    if draw == 'wire':
        model.Display(True)
    if draw == 'frame':
        model.Display(False)

    # Solve the model
    model.Analyze()

    # Find the deflection at the end
    N2 = model.GetNode('N2')
    delta_Z = N2.DZ * 1E3
    theoryF = (  force * length**3 / (3 * material.E * section.Iyy)) * 1E3
    theoryM = (-moment * length**2 / (2 * material.E * section.Iyy)) * 1E3
    theoryU = (uniform * length**4 / (8 * material.E * section.Iyy)) * 1E3
    print('cantilever: Deflection at end = {d:.3}mm (expected {e:.3g}mm from beam theory)'.format(d=delta_Z, e=(theoryF + theoryM + theoryU)))
    
    if plot:
        M1 = model.GetMember('M1')
        M1.PlotShear("Fz")
        M1.PlotMoment("My")
        M1.PlotDeflection("dz")

    if stress is not None: # Note: experimental
        model.DisplayResults(stress)

def three_point(section, material, length, force, moment, **kwargs):
    """Analyse a three-point bending problem
    """
    # Add some options to display frame & results
    if 'draw' in kwargs: # Note: requires VisPy [options: 'draw', 'frame']
        draw = kwargs['draw']
    else:
        draw = None
    if 'plot' in kwargs: # [options: True, False]
        plot = kwargs['plot']
    else:
        plot = False
    if 'stress' in kwargs: # Note: experimental; requires VisPy [options: 'seq', 'sxx', 'txy', 'tzx']
        stress = kwargs['stress']
    else:
        stress = None

    # Create a new finite element model                                                                                                                                                                                             
    model = FEModel3D()

    # Add nodes (N1, N2 and N3) at ends and middle of beam; offset from origin for the sake of visualisation                                                                                                                                        
    model.AddNode('N1', -length/2, 0.1, 0.1)
    model.AddNode('N2',         0, 0.1, 0.1)
    model.AddNode('N3',  length/2, 0.1, 0.1)

    # Add beams (M1 and M2) between the nodes
    model.AddMemberExt('M1', 'N1', 'N2', material, section)
    model.AddMemberExt('M2', 'N2', 'N3', material, section)

    # Simply supported - can't move, can't twist about X- or Z-axes
    model.DefineSupport('N1',  True, True, True, True, False, True)
    model.DefineSupport('N3', False, True, True, True, False, True) # Allow movement in X-direction

    if force != 0:
        # Add a vertical (Z) load at the end
        model.AddNodeLoad('N2', 'FZ', force)

    if moment != 0:
        # Add a moment at the end in the ZX-plane
        model.AddNodeLoad('N2', 'MY', moment)

    if draw == 'wire':
        model.Display(True)
    if draw == 'frame':
        model.Display(False)

    # Solve the model
    model.Analyze()

    # Find the deflection at the end
    N2 = model.GetNode('N2')
    delta_Z = N2.DZ * 1E3
    theoryF = (  force * length**3 / (48 * material.E * section.Iyy)) * 1E3
    theoryM = 0
    print('three_point: Deflection at mid-point = {d:.3}mm (expected {e:.3g}mm from beam theory)'.format(d=delta_Z, e=(theoryF + theoryM)))
    
    if plot:
        M1 = model.GetMember('M1')
        M1.PlotShear("Fz")
        M1.PlotMoment("My")
        M1.PlotDeflection("dz")

    if stress is not None: # Note: experimental
        model.DisplayResults(stress)

def bracket(section, material, length, force, **kwargs):
    """Analyse a three-member bracket as a pin-jointed truss
    """
    # Add some options to display frame & results
    if 'draw' in kwargs: # Note: requires VisPy [options: 'draw', 'frame']
        draw = kwargs['draw']
    else:
        draw = None
    if 'plot' in kwargs: # [options: True, False]
        plot = kwargs['plot']
    else:
        plot = False
    if 'stress' in kwargs: # Note: experimental; requires VisPy [options: 'seq', 'sxx', 'txy', 'tzx']
        stress = kwargs['stress']
    else:
        stress = None

    # Create a new finite element model                                                                                                                                                                                             
    model = FEModel3D()

    # Add nodes (N1, N2 and N3) at wall, with N4 as load point                                                                                                                                        
    model.AddNode('N1',      0, -length/2,  length/4)
    model.AddNode('N2',      0,  length/2,  length/2)
    model.AddNode('N3',      0,         0, -length/2)
    model.AddNode('N4', length,         0,  length/2)

    # Add beams (M1, M2 and M3) between the nodes
    model.AddMemberExt('M1', 'N1', 'N4', material, section)
    model.AddMemberExt('M2', 'N2', 'N4', material, section)
    model.AddMemberExt('M3', 'N3', 'N4', material, section)

    # Pin-jointed: Fix the support nodes...
    model.DefineSupport('N1', True, True, True, True, True, True)
    model.DefineSupport('N2', True, True, True, True, True, True)
    model.DefineSupport('N3', True, True, True, True, True, True)

    # but allow non-twist gradient at beam ends
    model.DefineReleases('M1', False, False, False, False, True, True,
                               False, False, False, False, True, True)
    model.DefineReleases('M2', False, False, False, False, True, True,
                               False, False, False, False, True, True)
    model.DefineReleases('M3', False, False, False, False, True, True,
                               False, False, False, False, True, True)

    if force != 0:
        # Add a vertical (Z) load at the end
        model.AddNodeLoad('N4', 'FZ', force)

    if draw == 'wire':
        model.Display(True)
    if draw == 'frame':
        model.Display(False)

    # Solve the model
    model.Analyze()

    # Find the deflection at the end
    N4 = model.GetNode('N4')
    delta_Z = N4.DZ * 1E3
    print('bracket: Deflection = {d:.3}mm'.format(d=delta_Z))

    if plot:
        M3 = model.GetMember('M3')
        M3.PlotShear("Fz")
        M3.PlotMoment("My")
        M3.PlotDeflection("dz")

    if stress is not None: # Note: experimental
        model.DisplayResults(stress)

#cantilever(section_6, material_2, 3, -1000,     0) # downwards force only
#cantilever(section_6, material_2, 3,     0, -1000) # clockwise moment only
#cantilever(section_6, material_2, 3,  -750, -2000, draw='wire')
#cantilever(section_6, material_2, 3, -1250, -2000, plot=True)
#cantilever(section_6, material_2, 3,     0,     0, gravity=True, plot=True)

#three_point(section_6, material_2, 3, -1000,     0, plot=True) # downwards force only
#three_point(section_6, material_2, 3,     0, -1000, plot=True) # clockwise moment only
#three_point(section_6, material_2, 3, -1000, -1000, plot=True) # clockwise moment only
#three_point(section_6, material_2, 3,     0, -1000, stress='sxx')

#bracket(section_4, material_1, 1, -10000, stress='seq')
