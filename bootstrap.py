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

		self.table = tbl[['ID','sfg','ra','dec','z_peak','z_err','LMASS','LMASS_ERR']]

	def perturb_catalog(self, perturb_z = False):

		pseudo_cat = self.table.copy()
		ngals = len(pseudo_cat)

		if perturb_z == True:
			pseudo_z = self.table['z_peak'] + self.table['z_err']*np.random.randn(len(self.table['z_err']))
			pseudo_cat['z_peak'] = pseudo_z

		self.pseudo_cat = pseudo_cat.sample(ngals,replace=True) 
		#pdb.set_trace()