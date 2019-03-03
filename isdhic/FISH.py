#### Visualisation of tads
# Example run: run /Users/agalicina/Google Drive/Lab_mac/SINGLE_CELL/MODELLING1/WD/script_probes.py
### Radius of hyration calculation


#python 

from pymol import cmd
import math

tad_command = """
select tmp1, ({cell} and id {bgn}-{end})
create tad_{name}, tmp1
"""

probe_command = """
select tmp1, ({cell} and id {pos})
create probe_{name}, tmp1
"""





cmd.do("""
load()

bg_color white

set transparency, 0.8
set solvent_radius, 1
set surface_quality, 1
set surface_color, gray
set sphere_scale, 0.1
set_color natred=[0.81, 0.24, 0.16]
set_color natorange=[0.95, 0.58, 0.19]
set_color natblue=[0.21, 0.61, 0.84]
set_color natpurple=[0.65, 0.36, 0.88]
set_color natgreen=[0.27, 0.60, 0.38]
""")



tads_coords = {
    'S1T1':[376, 386],
    'S1T2':[386, 400],
    'S1T3':[400, 411],
    'S1T0':[375, 412],
    
    'S2T1':[1760, 1764],
    'S2T2':[1764, 1772],
    'S2T3':[1772, 1778],
    'S2T0':[1759, 1779],
    
    'S3T1':[703, 718],
    'S3T2':[718, 722],
    'S3T3':[722, 740],
    'S3T4':[740, 755],
    'S3T0':[702, 756]
}

probes_coords = {
    'S1P1': 388,
    'S1P2': 397,
    'S1P3': 406,
    'S2P1':1765,
    'S2P2':1771,
    'S2P3':1777,
    'S3P1':724,
    'S3P2':736,
    'S3P3':714,
    'S3P4':747
}



experiment_ids = "A2".upper().split()

for name, (bgn, end) in tads_coords.iteritems():
	formatting_dct = {'cell': 'B19.sampled'}
	formatting_dct.update({'bgn':bgn+1, 'end':end, 'name':name})
        cmd.do(tad_command.format(**formatting_dct))

for name, pos in probes_coords.iteritems():
	formatting_d = {'cell': 'B19.sampled'}
	formatting_d.update({'pos':pos, 'name':name})
        cmd.do(probe_command.format(**formatting_d))


cmd.do("""
show cartoons, (tad*)
show spheres, (probe*)

set cartoon_color, gray80, (tad*T0)

set cartoon_color, natred, (tad*T1)
set cartoon_color, natblue, (tad*T2)
set cartoon_color, natorange, (tad*T3)
set cartoon_color, natgreen, (tad*T4)
set cartoon_color, natorange, (tad_S3T2)
set cartoon_color, natblue, (tad_S3T3)

color natred, (tad_*1)
color natblue, (tad_*2)
color natorange, (tad_*3)
color natgreen, (tad_*4)
color natorange, (tad_S3T2)
color natblue, (tad_S3T3)

color deepblue, (probe_*1)
color teal, (probe_*2)
color natorange, (probe_*3)
color natgreen, (probe_*4)
""")


