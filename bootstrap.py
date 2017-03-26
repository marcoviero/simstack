# Standard modules
import pdb
import os.path
import sys
import six
import shutil
import time
import logging
import importlib
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.wcs import WCS

class Bootstrap:

	def __init__(self,tbl):

		#pdb.set_trace()
		#try:
		#	self.table = tbl[['ID','sfg','ra','dec','z_peak','z_err','LMASS','LMASS_ERR','mips24','F_ratio']]
		#except:
		#	self.table = tbl[['ID','sfg','ra','dec','z_peak','LMASS','mips24','F_ratio']]
		self.table = tbl

	def perturb_catalog(self, perturb_z = False, draw_z = False):
		'''
		Default is a simple bootstrap with replacement.
		If perturb_z is True, then the redshifts (z_peak) are perturbed by random normal with sigma z_err
		If draw_z is True, then the redshifts are drawn from the redshift PDF (output by EAZY) of each object (TODO!)
		'''

		pseudo_cat = self.table.copy()
		ngals = len(pseudo_cat)

		if perturb_z == True:
			pseudo_z = self.table['z_peak'] + self.table['z_err']*np.random.randn(len(self.table['z_err']))
			pseudo_cat['z_peak'] = pseudo_z

		#Simple Bootstrap.  Sample draws from pseudo_cat with replacement.
		self.pseudo_cat = pseudo_cat.sample(ngals,replace=True)
		#pdb.set_trace()
