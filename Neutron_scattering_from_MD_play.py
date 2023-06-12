#!/usr/bin/env python
# coding: utf-8

import periodictable
from Bio.PDB import PDBParser
from refnx.reflect import SLD, Slab, Structure, ReflectModel
import MDAnalysis as mda
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np
import os.path

import refnx, scipy 
from refnx.reflect import reflectivity
# the ReflectDataset object will contain the data
from refnx.dataset import ReflectDataset
# the reflect module contains functionality relevant to reflectometry
from refnx.reflect import ReflectModel, SLD, Stack
# the analysis module contains the curvefitting engine
from refnx.analysis import Objective, Transform, CurveFitter

from md_simulation import MDSimulation


# See this for help: https://nbviewer.org/github/refnx/refnx-models/blob/master/md_simulation/md_simulation.ipynb 
# 
# 
# For some background knowledge: http://gisaxs.com/index.php/Scattering_Length_Density#:~:text=The%20Scattering%20Length%20Density%20(SLD,of%20the%20'scattering%20entities'
# 
# And here_: https://periodictable.readthedocs.io/en/latest/guide/intro.html

## A good reference as well : https://www.mdpi.com/2218-273X/12/11/1591


### for CG simulations refer to this thesis and paper:  https://www.proquest.com/openview/25f2890d3d96d6883c8c88310c9a2152/1?pq-origsite=gscholar&cbl=18750&diss=y
## https://pubs.acs.org/doi/full/10.1021/acs.jctc.0c00132



def strip_end(string):
    """
    Strips 'R', 'S', or 'T' from the end of a string if the string has more than 4 characters.
    
    Args:
        string (str): The input string.
        
    Returns:
        str: The stripped string.
    """
    if len(string) > 4:
        while string.endswith(('R', 'S', 'T')):
            string = string[:-1]
    
    return string



# #### Extracts the SLD values for each atom type in the universe



# Load trajectory and topology and strip for ions
u = mda.Universe("step7_0.gro", "step7_0.part0002.xtc")
sub_uni = u.select_atoms('all and not resname SOD CLA')

with mda.Writer("sub.xtc", sub_uni.n_atoms) as W:
    for ts in u.trajectory:
        W.write(sub_uni)
sub_uni.atoms.write('sub.pdb')

##Reload trajectory
u = mda.Universe('sub.pdb', 'sub.xtc')

# list of atoms to be kept as hydrogenous in the mc3 lipid 
mc3_H_atoms = 'HN11 HN12 HN13 HN21 HN22 HN23 H11 H12 H21 H22 H31 H32 H5\
                H11R H11S H12R H12S'

# for charged mc3
mc3_H_atoms = 'HN4 HN11 HN12 HN13 HN21 HN22 HN23 H11 H12 H21 H22 H31 H32 H5\
                H11R H11S H12R H12S'

keep_H = mc3_H_atoms + mc3_H_atoms

# for polyA
polyA_H_atoms = "H61 H62"

#for resname in np.unique(u.select_atoms('all').resnames):
for resname in ['DLMC', 'TIP3']:
	u = mda.Universe('sub.pdb', 'sub.xtc')
	if resname  == 'DLMC':
		keep_H =  mc3_H_atoms + mc3_H_atoms
		#print (f'resname {resname} and type H and not name {keep_H}')
		u.select_atoms(f'resname {resname} and type H and not name {keep_H}').atoms.types = 'D'
	if resname == 'ADE':
		keep_H = polyA_H_atoms
		#print (f'resname {resname} and type H and not name {keep_H}')
		u.select_atoms(f'resname {resname} and type H and not name {keep_H}').atoms.types = 'D'
	else:
		#print (f'resname {resname} and type H')
		u.select_atoms(f'resname {resname} and type H').atoms.types = 'D'
	
	
	sim = MDSimulation(u,flip=False, cut_off=0, layer_thickness=1, roughness=3.5)
	# this is just one of three ways to determine the scattering lengths
	sim.assign_scattering_lengths('neutron')
	sim.run()
	
	np.save(f'SLD_{resname}.npy', sim.sld_profile())
	plt.plot(sim.sld_profile()[0], sim.sld_profile()[1]*10, label=resname)
	plt.xlabel('$z$/Å')
	plt.ylabel(r'$\rho(z)$/Å$^{-2}$')


	"""
	generate reflectivity from MD-SLD with refnx 
	"""
	
	# changing to refnx notation for clarity 
	structure = sim.sld_profile()[1]
	
	# define roughness of each layer 
	roughness = 3.5
	thickness = 1.0
	
	# define imaginary SLD component 
	SLD_im = 0
	
	# define x (q) and y (R) variables 
	q = np.linspace(0.01,0.35,int(1e3))
	
	  
	# define empty array for populating with individual slab layers 
	slabs = []
	
	# generate slabs for each structure 
	for SLD_re in structure:
	    slabs.append([thickness, SLD_re, SLD_im, roughness])
	
	# convert to array for refnx processing 
	slabs = np.array(slabs)
	
	# calculate and store reflectivity 
	R = list( reflectivity(q,slabs,bkg=8.0e-6) )
	np.save(f'Reflectivity_{resname}.npy', R)
	#plt.plot(q,R)
	#plt.xlim(0,0.1)
	#plt.ylim(0,0.005)


plt.legend(loc='best')
plt.xlim(80,150)
plt.savefig('SLD_profile_AA.png', dpi=300, bbox_inches='tight')
plt.savefig('Ref_profile_AA.png', dpi=300, bbox_inches='tight')










