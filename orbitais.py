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

#	Manually enter this array of values, you should open the same xyz file that is 
#	loaded as coords with pymol. In a terminal (some linux distro):
#	pymol my_coords.xyz
#	Once in pymol, enter get_view on the command line e copy the result.
#	You should format your set of values as presented on "the_view" 
#	variable.

the_view = ( 
				0.994782925,   -0.010882307,    0.101430826,
				0.077204466,    0.730216444,   -0.678839207,
				-0.066679209,    0.683128655,    0.727247536,
				0.000000000,    0.000000000,  -19.934177399,
				-0.000000238,    0.000000238,   -0.000000333,
			-5344.827636719, 5384.695800781,  -20.000000000 
			)

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

def create_png(coords,cube=None,labelname=None,colorA='orange',colorB='tea'):
	"""
	This is a simple function to create nice png images from a xyz file
	and a set of cube files. If only the xyz file is passed, this 
	function will return a png image containing the molecular geometry 
	stored in it using PyMol's cmd. You should provide set_view, 
	open pymol and load the xyz file, type get_view at the prompt and 
	copy the values to "cmd.set_view". If the script find some cube 
	files on the work directory, it will try to render isosurface's 
	images. You may wish to change the isosurface's color, this can be 
	done by altering the values of colorA and colorB parameters. 
	"""

	cmd.load(coords,object='coords',format='xyz',quiet=1)
	cmd.set_view(the_view)
	
	#	It is possible to set isovalues, colors and other stuffs
	#	the values you see are just the ones I used to prefer
	if cube:	
		cmd.load(cube,'Orb',format='cube',state=0)
		cmd.isosurface('Asurf1','Orb', 0.03,state=0)
		cmd.isosurface('Bsurf1','Orb', -0.03,state=0)
		cmd.color(color=colorA,selection='Asurf1')
		cmd.set('transparency',0.0,'Asurf1')
		cmd.color(color=colorB,selection='Bsurf1')
		cmd.set('transparency',0.0,'Bsurf1')
		
		#	Label position is a headache, you should consider set the best 
		#	position for your visualization. GOOD LUCK!
		
		cmd.pseudoatom(object='foo',label=labelname,pos=[-3,3,1.75])
#		cmd.show('spheres','foo')
#		cmd.label('foo',labelname)
		cmd.set("label_size",-2)
		cmd.set("label_outline_color","black")
		cmd.set("label_color","black")
		cmd.set("label_font_id",7)
#		cmd.set("label_position",[-6,7,1.75])
		cmd.png(labelname,width=1980,height=1740,dpi=200,ray=1)
		cmd.delete('Orb') 
		if cube:
			cmd.delete('Orb')
			cmd.delete('Asurf1')
			cmd.delete('Bsurf1')
			cmd.delete('foo')
			cmd.png('coords',width=1980,height=1740,dpi=200,ray=1)	
	else:
		cmd.png('coords',width=1980,height=1740,dpi=200,ray=1)
	return 

if not file_om:
	print("Nothing to do about surfaces here!\nI can generate a nice geom for you instead.")
	setview = the_view
	create_png(coords,the_view)
	
else:
	print("Let's get some plots!\n But I'll give you a nice\n geom picture.")	
	for cube in file_om:
		temp = re.search('(?<=mo)[0-9]+',cube)
		labelname = temp.group(0)
		create_png(coords,cube,labelname)
		
