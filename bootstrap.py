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

class StackedFluxes:

	def __init__(self, config_path, config_file):

		self.path = config_path
		self.config_file = config_file
		self.params = self.get_parameters(config_path+config_file)

		indpop = np.argsort(np.array([i for i in self.params['populations'].values()]))
		self.pops = [self.params['populations'].keys()[i] for i in indpop]
		self.npops = len(self.pops)
		self.nz = len(self.params['bins']['z_nodes']) - 1
		self.nm = len(self.params['bins']['m_nodes']) - 1
		self.nw = len(self.params['map_files'])
		self.ind = np.argsort(np.array([self.params['wavelength'][wv] for wv in self.params['wavelength']]))
		self.maps = [self.params['wavelength'].keys()[i] for i in self.ind]
		self.wvs = [self.params['wavelength'].values()[i] for i in self.ind]
		self.z_nodes = self.params['bins']['z_nodes']
		self.m_nodes = self.params['bins']['m_nodes']
		z_m_keys = self.m_z_key_builder()
		self.z_keys = z_m_keys[0]
		self.m_keys = z_m_keys[1]
		#self.stack_keys = self.stack_key_builder()
		#self.bootstrap_flux_array = np.array()
		#self.bootstrap_error_array = {}
		#self.bootstrap_error_dict = {}
		#self.binsize_error_array  = {}

	def get_parameters(self, config_path):

		params = parameters.get_params(config_path)

		return params

	def read_fluxes(self):
		#if (config_file): config_file = shortname_boots + '.cfg'

		stacked_fluxes = np.zeros([self.nw,self.nz,self.nm,self.npops])
		slice_keys = self.slice_key_builder()

		for i in range(self.nz):
			zboot = slice_keys[i]
			z_suf = 'z_'+ self.z_keys[i]
			filename_stacks = 'simstack_flux_densities_'+ self.params['io']['shortname'] + '_' + zboot + '.p'
			if os.path.exists(self.path+filename_stacks):
				simstack = pickle.load( open( self.path + filename_stacks, "rb" ))
				for wv in range(self.nw):
					single_wv_stacks = simstack[1][zboot][self.maps[wv]]
					for j in range(self.nm):
						m_suf = 'm_' + self.m_keys[j]
						#for p_suf in self.pops:
						for p in range(self.npops):
							p_suf = self.pops[p]
							#z_0p054_0p176__m_10p0_10p5_sf
							key = clean_args(z_suf+'__'+ m_suf+ '_' + p_suf)
							stacked_fluxes[wv,i,j,p] = single_wv_stacks[key].value

		self.simstack_flux_array = stacked_fluxes

	def slice_key_builder(self):

		if self.params['bins']['bin_in_lookback_time']:
			z_nodes = self.params['bins']['t_nodes']
		else:
			z_nodes = self.params['bins']['z_nodes']
		nz = len(z_nodes) - 1

		return [str(z_nodes[i])+ '-' +str(z_nodes[i+1]) for i in range(nz)]
		#return ['simstack_flux_densities_'+ self.params['io']['shortname'] + '_' + str(z_nodes[i])+ '-' +str(z_nodes[i+1]) + '_boot_' for i in range(nz)]

	def m_z_key_builder(self):

		z_suf = []
		m_suf = []

		for i in range(self.nz):
			if self.params['bins']['bin_in_lookback_time']:
				z_suf.append('{:.3f}'.format(self.params['bins']['z_nodes'][i]) +'-'+ '{:.3f}'.format(self.params['bins']['z_nodes'][i+1]))
			else:
				z_suf.append(str(self.params['bins']['z_nodes'][i])+'-'+str(self.params['bins']['z_nodes'][i+1]))

		for j in range(self.nm):
			m_suf.append(str(self.params['bins']['m_nodes'][j])+'-'+str(self.params['bins']['m_nodes'][j+1]))

		return [z_suf, m_suf]

class ErrorBars:

	def __init__(self, config_path, config_file):

		self.path = config_path
		self.config_file = config_file
		self.params = self.get_parameters(config_path+config_file)

		indpop = np.argsort(np.array([i for i in self.params['populations'].values()]))
		self.pops = [self.params['populations'].keys()[i] for i in indpop]
		self.npops = len(self.pops)
		self.nboots = int(self.params['number_of_boots'])
		self.nz = len(self.params['bins']['z_nodes']) - 1
		self.nm = len(self.params['bins']['m_nodes']) - 1
		self.nw = len(self.params['map_files'])
		self.ind = np.argsort(np.array([self.params['wavelength'][wv] for wv in self.params['wavelength']]))
		self.maps = [self.params['wavelength'].keys()[i] for i in self.ind]
		self.wvs = [self.params['wavelength'].values()[i] for i in self.ind]
		self.z_nodes = self.params['bins']['z_nodes']
		self.m_nodes = self.params['bins']['m_nodes']
		z_m_keys = self.m_z_key_builder()
		self.z_keys = z_m_keys[0]
		self.m_keys = z_m_keys[1]
		self.slice_keys = self.slice_key_builder()
		#self.stack_keys = self.stack_key_builder()
		#self.bootstrap_flux_array = np.array()
		#self.bootstrap_error_array = {}
		#self.bootstrap_error_dict = {}
		#self.binsize_error_array  = {}

	def get_error_bar_array(self):
		print 'return a numpy array'

	def get_error_bar_dictionary(self):
		print 'return a dictionary with nodes for keys'

	def get_parameters(self, config_path):

		params = parameters.get_params(config_path)

		return params

	def read_bootstraps(self):
		#if (config_file): config_file = shortname_boots + '.cfg'

		bootstrap_fluxes = np.zeros([self.nw,self.nz,self.nm,self.npops,self.nboots])
		slice_keys = self.slice_key_builder()

		#for zboot in slice_keys:
		for i in range(self.nz):
			zboot = slice_keys[i]
			z_suf = 'z_'+ self.z_keys[i]
			for k in np.arange(self.nboots) + int(self.params['boot0']):
				filename_boots = 'simstack_flux_densities_'+ self.params['io']['shortname'] + '_' + zboot + '_boot_'+ str(k) + '.p'
				if os.path.exists(self.path+filename_boots):
					bootstack = pickle.load( open( self.path + filename_boots, "rb" ))
					#for wv in self.maps:
					for wv in range(self.nw):
						single_wv_stacks = bootstack[1][zboot][self.maps[wv]]
						for j in range(self.nm):
							m_suf = 'm_' + self.m_keys[j]
							#for p_suf in self.pops:
							for p in range(self.npops):
								p_suf = self.pops[p]
								#z_0p054_0p176__m_10p0_10p5_sf
								key = clean_args(z_suf+'__'+ m_suf+ '_' + p_suf)
								bootstrap_fluxes[wv,i,j,p,k] = single_wv_stacks[key].value

		self.bootstrap_flux_array = bootstrap_fluxes


	def is_bootstrap(self,config):
		return config['bootstrap']

	def slice_key_builder(self):

		if self.params['bins']['bin_in_lookback_time']:
			z_nodes = self.params['bins']['t_nodes']
		else:
			z_nodes = self.params['bins']['z_nodes']
		nz = len(z_nodes) - 1

		return [str(z_nodes[i])+ '-' +str(z_nodes[i+1]) for i in range(nz)]
		#return ['simstack_flux_densities_'+ self.params['io']['shortname'] + '_' + str(z_nodes[i])+ '-' +str(z_nodes[i+1]) + '_boot_' for i in range(nz)]

	def m_z_key_builder(self):

		z_suf = []
		m_suf = []

		for i in range(self.nz):
			if self.params['bins']['bin_in_lookback_time']:
				z_suf.append('{:.3f}'.format(self.params['bins']['z_nodes'][i]) +'-'+ '{:.3f}'.format(self.params['bins']['z_nodes'][i+1]))
			else:
				z_suf.append(str(self.params['bins']['z_nodes'][i])+'-'+str(self.params['bins']['z_nodes'][i+1]))

		for j in range(self.nm):
			m_suf.append(str(self.params['bins']['m_nodes'][j])+'-'+str(self.params['bins']['m_nodes'][j+1]))

		return [z_suf, m_suf]


class PickledStacksReader:

	def __init__(self, config_path, config_file):
		''' Uses the config_file to determine if it is bootstraps or not'''

		self.path = config_path
		self.config_file = config_file
		self.params = self.get_parameters(config_path+config_file)
		if self.params['bootstraps']:
			self.nboots = int(self.params['number_of_boots'])

		indpop = np.argsort(np.array([i for i in self.params['populations'].values()]))
		self.pops = [self.params['populations'].keys()[i] for i in indpop]
		self.npops = len(self.pops)
		self.nz = len(self.params['bins']['z_nodes']) - 1
		self.nm = len(self.params['bins']['m_nodes']) - 1
		self.nw = len(self.params['map_files'])
		self.ind = np.argsort(np.array([self.params['wavelength'][wv] for wv in self.params['wavelength']]))
		self.maps = [self.params['wavelength'].keys()[i] for i in self.ind]
		self.wvs = [self.params['wavelength'].values()[i] for i in self.ind]
		self.z_nodes = self.params['bins']['z_nodes']
		self.m_nodes = self.params['bins']['m_nodes']
		z_m_keys = self.m_z_key_builder()
		self.z_keys = z_m_keys[0]
		self.m_keys = z_m_keys[1]
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

		if self.params['boostrap'] == True:
			bootstrap_fluxes = np.zeros([self.nw,self.nz,self.nm,self.npops,self.nboots])
		else:
			stacked_fluxes = np.zeros([self.nw,self.nz,self.nm,self.npops])
		slice_keys = self.slice_key_builder()

		#for zboot in slice_keys:
		for i in range(self.nz):
			z_slice = slice_keys[i]
			z_suf = 'z_'+ self.z_keys[i]
			for wv in range(self.nw):
				for j in range(self.nm):
					m_suf = 'm_' + self.m_keys[j]
					if self.params['boostrap'] == True:
						for k in np.arange(self.nboots) + int(self.params['boot0']):
							filename_boots = 'simstack_flux_densities_'+ self.params['io']['shortname'] + '_' + z_slice + '_boot_'+ str(k) + '.p'
							if os.path.exists(self.path+filename_boots):
								bootstack = pickle.load( open( self.path + filename_boots, "rb" ))
								single_wv_stacks = bootstack[1][z_slice][self.maps[wv]]
								for p in range(self.npops):
									p_suf = self.pops[p]
									key = clean_args(z_suf+'__'+ m_suf+ '_' + p_suf)
									bootstrap_fluxes[wv,i,j,p,k] = single_wv_stacks[key].value
					else:
						filename_stacks = 'simstack_flux_densities_'+ self.params['io']['shortname'] + '_' + z_slice + '.p'
						if os.path.exists(self.path+filename_stacks):
							simstack = pickle.load( open( self.path + filename_stacks, "rb" ))
							single_wv_stacks = simstack[1][z_slice][self.maps[wv]]
							for p in range(self.npops):
								p_suf = self.pops[p]
								key = clean_args(z_suf+'__'+ m_suf+ '_' + p_suf)
								stacked_fluxes[wv,i,j,p,k] = single_wv_stacks[key].value

		if self.params['boostrap'] == True:
			self.bootstrap_flux_array = bootstrap_fluxes
		else:
			self.stacked_flux_array = stacked_fluxes

	def is_bootstrap(self,config):
		return config['bootstrap']

	def slice_key_builder(self):

		if self.params['bins']['bin_in_lookback_time']:
			z_nodes = self.params['bins']['t_nodes']
		else:
			z_nodes = self.params['bins']['z_nodes']
		nz = len(z_nodes) - 1

		return [str(z_nodes[i])+ '-' +str(z_nodes[i+1]) for i in range(nz)]
		#return ['simstack_flux_densities_'+ self.params['io']['shortname'] + '_' + str(z_nodes[i])+ '-' +str(z_nodes[i+1]) + '_boot_' for i in range(nz)]

	def m_z_key_builder(self):

		z_suf = []
		m_suf = []

		for i in range(self.nz):
			if self.params['bins']['bin_in_lookback_time']:
				z_suf.append('{:.3f}'.format(self.params['bins']['z_nodes'][i]) +'-'+ '{:.3f}'.format(self.params['bins']['z_nodes'][i+1]))
			else:
				z_suf.append(str(self.params['bins']['z_nodes'][i])+'-'+str(self.params['bins']['z_nodes'][i+1]))

		for j in range(self.nm):
			m_suf.append(str(self.params['bins']['m_nodes'][j])+'-'+str(self.params['bins']['m_nodes'][j+1]))

		return [z_suf, m_suf]
