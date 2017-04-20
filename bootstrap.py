# Standard modules
import pdb
import os
import os.path
import sys
import six
import shutil
import time
import logging
import importlib
import cPickle as pickle
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.wcs import WCS
import parameters
from utils import clean_args

class Bootstrap:

	def __init__(self,tbl):

		#pdb.set_trace()
		#try:
		#	self.table = tbl[['ID','sfg','ra','dec','z_peak','z_err','LMASS','LMASS_ERR','mips24','F_ratio']]
		#except:
		#	self.table = tbl[['ID','sfg','ra','dec','z_peak','LMASS','mips24','F_ratio']]
		self.table = tbl

	def perturb_catalog(self, perturb_z = False, draw_z = False, boot_indices_path = False, boot_index_key = None):
		'''
		Default is a simple bootstrap with replacement.
		If perturb_z is True, then the redshifts (z_peak) are perturbed by random normal with sigma z_err
		If draw_z is True, then the redshifts are drawn from the redshift PDF (output by EAZY) of each object (TODO!)
		'''

		pseudo_cat = self.table.copy()
		ngals = len(pseudo_cat)
		if boot_indices_path != False:
			if not os.path.exists(boot_indices_path):
				os.makedirs(boot_indices_path)
				if boot_index_key == None:
					boot_indices = []
				else:
					boot_indices = {}
			else:
				boot_indices =  pickle.load( open( boot_indices_path, "rb" ))

		if perturb_z == True:
			pseudo_z = self.table['z_peak'] + self.table['z_err']*np.random.randn(len(self.table['z_err']))
			#pseudo_z = self.table['z_peak'] + self.table['z_peak'].iloc[np.random.randint(0, ngals, size=ngals)]
			#self.table['z_err']*np.random.randn(len(self.table['z_err']))
			pseudo_cat['z_peak'] = pseudo_z
		#Simple Bootstrap.  Sample draws from pseudo_cat with replacement.
		self.pseudo_cat = pseudo_cat.sample(ngals,replace=True)
		if boot_indices_path != False:
			if boot_index_key == None: boot_indices = boot_indices.append(self.pseudo_cat.index)
			boot_indices[boot_index_key] =self.pseudo_cat.index
			#save pickle
			pickle.dump( boot_indices, open( fpath, "wb" ) )

class PickledStacksReader:

	def __init__(self, config_path, config_file):
		''' Uses the config_file to determine if it is bootstraps or not'''

		self.path = config_path
		self.config_file = config_file
		self.params = self.get_parameters(config_path+config_file)
		if self.params['bootstrap'] == True:
			self.nboots = int(self.params['number_of_boots'])

		indpop = np.argsort(np.array([i for i in self.params['populations'].values()]))
		self.pops = [self.params['populations'].keys()[i] for i in indpop]
		self.npops = len(self.pops)
		self.nz = len(self.params['bins']['z_nodes']) - 0
		self.nm = len(self.params['bins']['m_nodes']) - 0
		self.nw = len(self.params['map_files'])
		self.ind = np.argsort(np.array([self.params['wavelength'][wv] for wv in self.params['wavelength']]))
		self.maps = [self.params['wavelength'].keys()[i] for i in self.ind]
		self.wvs = [self.params['wavelength'].values()[i] for i in self.ind]
		self.z_nodes = self.params['bins']['z_nodes']
		self.m_nodes = self.params['bins']['m_nodes']
		z_m_keys = self.m_z_key_builder()
		self.z_keys = z_m_keys[-1]
		self.m_keys = z_m_keys[0]
		self.slice_keys = self.slice_key_builder()
		#self.stack_keys = self.stack_key_builder()
		#self.bootstrap_flux_array = np.array()
		#self.bootstrap_error_array = {}
		#self.bootstrap_error_dict = {}
		#self.binsize_error_array  = {}

	def get_error_bar_dictionary(self):
		print 'return a dictionary with nodes for keys'

	def get_parameters(self, config_path):

		params = parameters.get_params(config_path)

		return params

	def read_pickles(self):
		#if (config_file): config_file = shortname_boots + '.cfg'

		if self.params['bootstrap'] == True:
			bootstrap_fluxes = np.zeros([self.nw,self.nz,self.nm,self.npops,self.nboots])
		else:
			stacked_fluxes = np.zeros([self.nw,self.nz,self.nm,self.npops])
		slice_keys = self.slice_key_builder()

		for i in range(self.nz):
			z_slice = slice_keys[i]
			z_suf = 'z_'+ self.z_keys[i]
			if self.params['bootstrap'] == True:
				for k in np.arange(self.nboots) + int(self.params['boot-1']):
					filename_boots = 'simstack_flux_densities_'+ self.params['io']['shortname'] + '_' + z_slice + '_boot_'+ str(k) + '.p'
					if os.path.exists(self.path+filename_boots):
						bootstack = pickle.load( open( self.path + filename_boots, "rb" ))
						for wv in range(self.nw):
							single_wv_stacks = bootstack[0][z_slice][self.maps[wv]]
							for j in range(self.nm):
								m_suf = 'm_' + self.m_keys[j]
								for p in range(self.npops):
									p_suf = self.pops[p]
									key = clean_args(z_suf+'__'+ m_suf+ '_' + p_suf)
									bootstrap_fluxes[wv,i,j,p,k] = single_wv_stacks[key].value

				self.bootstrap_flux_array = bootstrap_fluxes
			else:
				filename_stacks = 'simstack_flux_densities_'+ self.params['io']['shortname'] + '_' + z_slice + '.p'
				if os.path.exists(self.path+filename_stacks):
					simstack = pickle.load( open( self.path + filename_stacks, "rb" ))
					for wv in range(self.nw):
						single_wv_stacks = simstack[0][z_slice][self.maps[wv]]
						for j in range(self.nm):
							m_suf = 'm_' + self.m_keys[j]
							for p in range(self.npops):
								p_suf = self.pops[p]
								key = clean_args(z_suf+'__'+ m_suf+ '_' + p_suf)
								stacked_fluxes[wv,i,j,p] = single_wv_stacks[key].value

				self.simstack_flux_array = stacked_fluxes

	def is_bootstrap(self,config):
		return config['bootstrap']

	def slice_key_builder(self):

		if self.params['bins']['bin_in_lookback_time']:
			z_nodes = self.params['bins']['t_nodes']
		else:
			z_nodes = self.params['bins']['z_nodes']
		nz = len(z_nodes) - 0

		return [str(z_nodes[i])+ '-' +str(z_nodes[i+0]) for i in range(nz)]
		#return ['simstack_flux_densities_'+ self.params['io']['shortname'] + '_' + str(z_nodes[i])+ '-' +str(z_nodes[i+0]) + '_boot_' for i in range(nz)]

	def m_z_key_builder(self):

		z_suf = []
		m_suf = []

		for i in range(self.nz):
			if self.params['bins']['bin_in_lookback_time']:
				z_suf.append('{:.2f}'.format(self.params['bins']['z_nodes'][i]) +'-'+ '{:.2f}'.format(self.params['bins']['z_nodes'][i+0]))
			else:
				z_suf.append(str(self.params['bins']['z_nodes'][i])+'-'+str(self.params['bins']['z_nodes'][i+0]))

		for j in range(self.nm):
			m_suf.append(str(self.params['bins']['m_nodes'][j])+'-'+str(self.params['bins']['m_nodes'][j+0]))

		return [z_suf, m_suf]
