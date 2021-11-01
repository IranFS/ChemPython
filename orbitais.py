#!/usr/bin/env python3.9
# -*- coding: utf-8 -*-
##################################---Description----######################################
#	A very rudimentar and initial iteration of a script that stores the cubes inside a 
# folder and, given that the corresponding xyz file is entered as an argument, generates
# quite beautiful png files for the surfaces of each cube file.
# 
#	It is a, at some extent, hardcoded script. If you find it and knows how to improve
# it, fell free. Also, share your improvements. 
#
# Partially based on https://gist.github.com/bobbypaton/1cdc4784f3fc8374467bae5eb410edef
# thanks, dude!
##########################################################################################

from glob import glob
import os
import sys
import re

# setting pymol's environment

moddir='/home/softwares/pymol/lib/python/'
sys.path.insert(0, moddir)
os.environ['/home/softwares/pymol/bin'] = os.path.join(moddir, 'pymol/pymol_path')

# import pymol and cmd

import pymol
pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
pymol.finish_launching()
cmd = pymol.cmd

# acquire the xyz file and the cube files inside the folder
# note that the xyz file must be passed as an argument

coords = sys.argv[1]
file_om=[f for f in glob("*.cube")]
print(file_om)
Oms = len(file_om)

# workspace settings

cmd.bg_color("white")
cmd.set("ray_opaque_background", "off")
cmd.set("field_of_view", 100)
cmd.set("orthoscopic", 'off')
cmd.set("transparency", 0)
cmd.set("dash_gap", 0)
cmd.set("ray_trace_mode", 1)
cmd.set("ray_trace_fog",0)
cmd.set("ray_texture", 2)
cmd.set("antialias", 1)
cmd.set("ambient", 1)
cmd.set("spec_count", 5)
cmd.set("shininess", 90)
cmd.set("specular", 1)
cmd.set("reflect", .1)
cmd.space("cmyk")
cmd.show("sticks")
cmd.show("spheres")
cmd.color("gray85","elem C")
cmd.color("gray98","elem H")
cmd.color("slate","elem N")
cmd.set("stick_radius",0.07)
cmd.set("sphere_scale",0.18)
cmd.set("sphere_scale",0.13, "elem H")
cmd.set("dash_gap",0.01,)
cmd.set("dash_radius",0.07,)
cmd.set("stick_color","black",)
cmd.set("dash_gap",0.01)
cmd.set("dash_radius",0.035)

# Set to build bonds between Europium and coordinated atoms on its complexes
# if the desired element is not present, include a similar line, set the 
# bond lenght also

cmd.distance('coordinate','elem Eu', 'elem N', 2.7,label=0)
cmd.distance('coordinate','elem Eu', 'elem O', 2.7,label=0)


# Start a rudimentary for loop that iterates over all the files 
# inside file_om and get the regex inside .mo****a. this is a far 
# more soft solution than try to slice each string in file_om 
# with the file_om[i][:] method 

for cube in file_om:
	temp = re.search('(?<=mo)[0-9]+',cube)
	labelname = temp.group(0)

#	This section performs a small iteration over the list that stores the .cube files 
#	parsed by glob, for each cube file a sequency of PyMol's cmd calls are made, and 
#	the loop will end generationg a nice png file

for cube in file_om:
	temp = re.search('(?<=mo)[0-9]+',cube)
	labelname = temp.group(0)
	cmd.load(coords,object='coords',format='xyz',quiet=1)

#	Manuaally enter this array of values, you should open the same xyz file that is 
#	loaded as coords with pymol. In a terminal (some linux distro):
#	pymol my_coords.xyz
#	Once in pymol, enter get_view on the command line e copy the result.
	cmd.set_view((\
				  0.414068997,    0.909954250,    0.023042712,\
				  0.329461694,   -0.173420861,    0.928105175,\
				  0.848530829,   -0.376706958,   -0.371601909,\
				  0.000000000,    0.000000000,  -43.606266022,\
				  1.609671593,   -0.394682407,   -0.436405420,\
				  39.733062744,   47.479469299,  -20.000000000 ))
  
#	It is possible to set isovalues, colors and other stuffs
#	the values you see are just the ones I used to prefer

	cmd.load(cube,'Orb',format='cube',state=0)
	cmd.isosurface('Asurf1','Orb', 0.03,state=0)
	cmd.isosurface('Bsurf1','Orb', -0.03,state=0)
	cmd.color(color="yellow",selection='Asurf1')
	cmd.set('transparency',0.0,'Asurf1')
	cmd.color(color="blue",selection='Bsurf1')
	cmd.set('transparency',0.0,'Bsurf1')
	
#	Label position is a headache, you should consider set the best 
#	position for your visualization. GOOD LUCK!

	cmd.pseudoatom(object='foo',pos=[0,-9,1.75])
	cmd.show('spheres','foo')
	cmd.label('foo',labelname)
	cmd.set("label_size",-2)
	cmd.set("label_outline_color","black")
	cmd.set("label_color","black")
	cmd.set("label_font_id",10)
	cmd.set("label_position",[-3,4,1.75])
	cmd.png(labelname,width=660,height=580,ray=1)
	cmd.delete('Orb') 

