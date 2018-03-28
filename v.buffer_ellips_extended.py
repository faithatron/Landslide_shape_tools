#!/usr/bin/env python

############################################################################
#
# MODULE:	   v.buffer_ellips.py
#
# AUTHOR(S):	Faith Taylor		   
# PURPOSE:	  Simulate road blockages due to random landslide triggered events
#
# DATE:		 Mon Oct 27 11:34:20 2014
# COPYRIGHT:	(C) 2014  Faith Taylor
#
#	   This program is free software under the GNU General Public
#	   License (>=v2). Read the file COPYING that comes with GRASS
#	   for details.
#
#
#############################################################################

#%module
#% description: Create elliptical lansldide shapes 
#% keywords: vector
#% keywords: landslide
#%end
#%option
#% key: true_ls
#% type: string
#% gisprompt: old,vector,vector
#% description: The real landslide inventory
#% required: yes
#%end
#%option
#% key: inventory_type
#% type: string
#% description: The type of shapes we are fitting to the landslides
#% required: yes
#%end
#%option
#% key: filepath
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% label: Name of input file to be imported
#% description: '-' for standard input
#% gisprompt: old,file,file
#%end

# Import necessary modules
from __future__ import division
import time
import math
from numpy import *
import time
import shlex
import sys
import os
import atexit

import grass.script as grass

from grass.script import core as grasscore


def cleanup():
	pass

def main():
	true_ls = options['true_ls']
	inventory_type = options['inventory_type']
	filepath = options['filepath']

	
	print inventory_type

	#create column headers
	columns = ['orig_ls_id integer',
			'east double precision',
			'north double precision',
			'radiusL double precision',
			'radiusW double precision',
			'azimuth double precision',
			'area_ls double precision',
			'perim_ls double precision']
			
	

# Create a vector map of randomly generated points. The accompanying table contains a landslide area and length to width ratio for each point. 
	grass.run_command('v.in.ascii', cat = 0, x = 2, y = 3, separator = 'comma',
			#input = "/home/faithtaylor/Desktop/Guat_mitch_NAD83_quad_elip_of_l_for_mapping.csv", 
			input = filepath,
			output = "real_ls_centroids",
			columns = columns,
			skip = 1,
			overwrite = True)
				
	grass.run_command('v.db.addcolumn', 
			map = "real_ls_centroids",
			columns = "area_elip double precision, area_elip_clip double precision, ellipticity double precision")
								  
	category_vals = grass.read_command("v.category",
			input = "real_ls_centroids",
			type = "point",
			option = "print",
			flags = "g")
	
							  						  
	category_vals_a = category_vals.split('\n')
			
	num_points = len(category_vals_a)
					  
	patch_list = []
							  
	for i in range(0, num_points-1):


		category_value = category_vals_a[i]
				
		#Extract the individual point, which will become the landslide centroid
		grass.run_command("v.extract", 
			input = "real_ls_centroids",
			output = "point"+str(category_value),
			cats = category_value,
			overwrite = True)
			
									
		#Read the major & minor axes and the rotation (aspect) of the landslide from the point file	
		major = grass.read_command("v.db.select",
					  				  flags = "c",
					  				  map = "point"+str(category_value), 
					  				  columns = "radiusL")
					  				  
		minor = grass.read_command("v.db.select",
					  				  flags = "c",
					  				  map = "point"+str(category_value), 
					  				  columns = "radiusW")	
	  				  						  				  		  
		aspect = grass.read_command("v.db.select",
					  				  flags = "c",
					  				  map = "point"+str(category_value), 
					  				  columns = "azimuth")

		#Add a buffer around that point to give the landslide an area
		grass.run_command("v.buffer",
						overwrite = True,
						input = "point"+str(category_value),
						output = str(inventory_type) +"Ellipse_landslide" +str(category_value),
						distance = major,
						minordistance = minor,
						angle = aspect,
						flags = "t")						

		# Verify the area of this ellipse
		grass.run_command("v.to.db",
				map = str(inventory_type) +"Ellipse_landslide" +str(category_value),
				option = "area",
				columns = "area_elip") 

			
		#Extract the matching real landslide
		grass.run_command("v.extract",
				input = true_ls,
				output = "Real_landslide"+str(category_value),
				cats = str(category_value),
				overwrite = True)		
		
		#Overlay this elliptical approximation of the landslide with the original landslide shape, where the category values should match
		grass.run_command("v.overlay", 
						overwrite = True,
						ainput = "Real_landslide"+str(category_value),
						binput = str(inventory_type) +"Ellipse_landslide"+str(category_value),
						operator = "and",
						output = str(inventory_type) +"AE_clip"+str(category_value))
						
		# Upload the area of this overlap to the file
		grass.run_command("v.to.db",
				map = str(inventory_type) +"AE_clip"+str(category_value),
				option = "area",
				columns = "b_area_elip_clip") 
				
		#Calculate the ellipticity of original shape
		grass.run_command("v.db.update",
			map = str(inventory_type) +"AE_clip"+str(category_value),
			column = "b_ellipticity",
			value = "1-((2*(b_area_ls - b_area_elip_clip))/b_area_ls)")
			
	patch_list = []
	for f in range (1, num_points):
		patch_name = (str(inventory_type)+'Ellipse_landslide'+str(f)+",")
		patch_list.append(patch_name)


	len_patch_list = len(patch_list)
	num_patches = int((len_patch_list/100))

	
	super_patch_list = []
	
	for i in range(0, num_patches):
		patch_list_a  = patch_list[(i*100):((i*100)+100)]
		patch_list_2 = ''.join(patch_list_a)
	
		if len(patch_list_a) > 1:
			
			
			grass.run_command("v.patch",
							input = ( patch_list_2),
							output = "inventory_patch"+str(i),
							flags = "e",
							overwrite = True)
	 
			super_patch_list.append("inventory_patch"+str(i)+",")
		else:
			super_patch_list.append(patch_list_2)
	
	#super_patch_list_2 = ''.join(super_patch_list)
	super_patch_list_2 = []
	
	for k in range(1, len(super_patch_list)):
		ind_patch_name = ("inventory_patch"+str(k)+",")
		super_patch_list_2.append(ind_patch_name)
	
	grass.run_command("v.patch",
				input = ( super_patch_list_2),
				output = str(inventory_type) +"Ellipse_inventory_final",
				flags = "e",
				overwrite = True)
####				
	patch_list = []
	for f in range (1, num_points):
		patch_name = (str(inventory_type)+'AE_clip'+str(f)+",")
		patch_list.append(patch_name)


	len_patch_list = len(patch_list)
	num_patches = int((len_patch_list/100))

	
	super_patch_list = []
	
	for i in range(0, num_patches+1):
		patch_list_a  = patch_list[(i*100):((i*100)+100)]
		patch_list_2 = ''.join(patch_list_a)
	
		if len(patch_list_a) > 1:
			
			
			grass.run_command("v.patch",
							input = ( patch_list_2),
							output = "clip_inventory_patch"+str(i),
							flags = "e",
							overwrite = True)
	 
			super_patch_list.append("clip_inventory_patch"+str(i)+",")
		else:
			super_patch_list.append(patch_list_2)
	
	#super_patch_list_2 = ''.join(super_patch_list)
	super_patch_list_2 = []
	
	for k in range(1, int((math.ceil((num_patches/100))))):
		ind_patch_name = ("clip_inventory_patch"+str(k)+",")
		super_patch_list_2.append(ind_patch_name)
	
	grass.run_command("v.patch",
				input = ( super_patch_list_2),
				output = str(inventory_type)+"Clip_inventory_final",
				flags = "e",
				overwrite = True)

if __name__ == "__main__":
	options, flags = grass.parser()
	atexit.register(cleanup)
	sys.exit(main())
			
				
