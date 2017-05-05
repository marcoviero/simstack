import pdb
import numpy as np
import gc
import os
import os.path
import sys
import cPickle as pickle
from astropy.wcs import WCS
import parameters
from utils import circle_mask
from utils import clean_args
from utils import clean_nans
from utils import dist_idl
from utils import gauss_kern
from utils import pad_and_smooth_psf
from utils import shift_twod
from utils import smooth_psf
from lmfit import Parameters, minimize, fit_report

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
		self.bin_ids = {}
		self.read_pickles()
		if self.params['bootstrap'] == True:
			ax = len(np.shape(self.bootstrap_flux_array)) - 1
			self.boot_error_bars = np.sqrt(np.var(self.bootstrap_flux_array,axis=ax))

			#self.covariance =

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
				for k in np.arange(self.nboots) + int(self.params['boot0']):
					filename_boots = 'simstack_flux_densities_'+ self.params['io']['shortname'] + '_' + z_slice + '_boot_'+ str(k) + '.p'
					if os.path.exists(self.path+filename_boots):
						bootstack = pickle.load( open( self.path + filename_boots, "rb" ))
						for bbk in bootstack[0]:
							self.bin_ids[bbk+'_'+str(k)] = bootstack[0][bbk]
						for wv in range(self.nw):
							try:
								single_wv_stacks = bootstack[1][z_slice][self.maps[wv]]
							except:
								single_wv_stacks = bootstack[1][self.maps[wv]]
							for j in range(self.nm):
								m_suf = 'm_' + self.m_keys[j]
								for p in range(self.npops):
									p_suf = self.pops[p]
									key = clean_args(z_suf+'__'+ m_suf+ '_' + p_suf)
									bootstrap_fluxes[wv,i,j,p,k] = single_wv_stacks[key].value
									#if wv ==0:
									#	self.bin_ids[key+'_'+str(k)] = bootstack[0]

				self.bootstrap_flux_array = bootstrap_fluxes
			else:
				filename_stacks = 'simstack_flux_densities_'+ self.params['io']['shortname'] + '_' + z_slice + '.p'
				if os.path.exists(self.path+filename_stacks):
					simstack = pickle.load( open( self.path + filename_stacks, "rb" ))
					#pdb.set_trace()
					#self.bin_ids = simstack[0]
					for ssk in simstack[0]:
						self.bin_ids[ssk] = simstack[0][ssk]
					for wv in range(self.nw):
						try:
							single_wv_stacks = simstack[1][z_slice][self.maps[wv]]
						except:
							single_wv_stacks = simstack[1][self.maps[wv]]
						for j in range(self.nm):
							m_suf = 'm_' + self.m_keys[j]
							for p in range(self.npops):
								p_suf = self.pops[p]
								key = clean_args(z_suf+'__'+ m_suf+ '_' + p_suf)
								stacked_fluxes[wv,i,j,p] = single_wv_stacks[key].value
								#if wv ==0:
								#	self.bin_ids[key] = simstack[0]
								#	pdb.set_trace()

				self.simstack_flux_array = stacked_fluxes

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


def simultaneous_stack_array_oned(p, layers_1d, data1d, err1d = None, arg_order = None):
  ''' Function to Minimize written specifically for lmfit '''

  v = p.valuesdict()

  len_model = len(data1d)
  nlayers = len(layers_1d)/len_model

  model = np.zeros(len_model)

  for i in range(nlayers):
    #print v.keys()[i]
    #print arg_order[i]
    if arg_order != None:
      model[:] += layers_1d[i*len_model:(i+1)*len_model] * v[arg_order[i]]
    else:
      model[:] += layers_1d[i*len_model:(i+1)*len_model] * v[v.keys()[i]]

  if err1d is None:
    return (data1d - model)
  return (data1d - model)/err1d

def simultaneous_stack_multimap(p, layers_1d, data1d, err1d = None, nmaps = 1, lmaps = None, arg_order = None):
  ''' Function to Minimize written specifically for lmfit '''

  flux = p.valuesdict()

  model = np.zeros(len(data1d))
  nlayers = len(layers_1d)/len(data1d)
  #if lmaps != None:
  #  imaps = np.cumsum([0] + lmaps)
  #else:
  #  imaps = [0]
  if nmaps == 1:
    imaps = [0,len(layers_1d)*nlayers]
  else:
    imaps = np.cumsum([0] + lmaps)

  for i in range(nmaps):
    i0 = imaps[i]*nlayers

    if nmaps == 1:
      len_model = len(data1d)
    else:
      len_model = lmaps[i]

    for j in range(nlayers):
      #print flux.keys()[i]
      #pdb.set_trace()
      if arg_order != None:
        model[imaps[i]:imaps[i+1]] += layers_1d[i0 + j*len_model:i0 + (j+1)*len_model] * flux[arg_order[j]]
      else:
        model[imaps[i]:imaps[i+1]] += layers_1d[i0 + j*len_model:i0 + (j+1)*len_model] * flux[flux.keys()[j]]

  if err1d is None:
    return (data1d - model)
  return (data1d - model)/err1d

def simultaneous_stack_array(p, layers_2d, data, err = None):
  ''' Function to Minimize written specifically for lmfit '''

  flux = p.valuesdict()

  csize = np.shape(layers_2d)

  model = np.zeros(csize[1])

  for i in range(csize[0]):
    #model += layers_2d[i,:] * flux['layer'+str(i)]
    model += layers_2d[i,:] * flux[flux.keys()[i]]

  if err is None:
    return (data - model)
  return (data - model)/err

def stack_in_redshift_slices(
  cmap,
  hd,
  layers_radec,
  fwhm=None,
  psf_names=None,
  cnoise=None,
  mask=None,
  beam_area=None,
  err_ss=None,
  quiet=None):
  ''' The first iteration of the translation from IDL to Python.
      Looks like an IDL function.
      Suggest using wrappers like viero_quick_stack.py
      but highly recommend Pythonic: stack_libraries_in_redshift_slices
      function that can be found below.

      Inputs:

  '''

  w = WCS(hd)
  #FIND SIZES OF MAP AND LISTS
  cms = np.shape(cmap)
  zeromask = np.zeros(cms)

  size_cube = np.shape(layers_radec)
  nsrcmax = size_cube[0]
  nlists = int(size_cube[1])

  ind_map_zero = np.where(np.isnan(cmap))
  nzero = np.shape(ind_map_zero)[1]

  if np.sum(cnoise) == 0: cnoise=cmap*0.0 + 1.0

  pix=hd["CD2_2"]*3600.
  if pix == 0: pix=hd["CDELT2"]*3600.

  #[STEP 0] - Calibrate maps
  if beam_area != None:
    cmap=cmap*beam_area*1e6
    cnoise=noise*beam_area*1e6

  # STEP 1  - Make Layers Cube
  layers=np.zeros([nlists,cms[0],cms[1]])

  for s in range(nlists):
    ind_src = np.where(layers_radec[:,s,0] != 0)
    if np.shape(ind_src)[1] > 0:
      ra = layers_radec[ind_src,s,0]
      dec = layers_radec[ind_src,s,1]
      # CONVERT FROM RA/DEC to X/Y
      # DANGER!!  NOTICE THAT I FLIP X AND Y HERE!!
      ty,tx = w.wcs_world2pix(ra, dec, 0)
      # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
      ind_keep = np.where((tx[0] >= 0) & (np.round(tx[0]) < cms[0]) & (ty[0] >= 0) & (np.round(ty[0]) < cms[1]))
      nt0 = np.shape(ind_keep)[1]
      real_x=np.round(tx[0,ind_keep][0]).astype(int)
      real_y=np.round(ty[0,ind_keep][0]).astype(int)
      # CHECK FOR SOURCES THAT FALL ON ZEROS MAP
      if nzero > 0:
        tally = np.zeros(nt0)
        for d in range(nt0):
          if cmap[real_x[d],real_y[d]] != 0:
            tally[d]=1.
        ind_nz=np.where(tally == 1)
        nt = np.shape(ind_nz)[1]
        real_x = real_x[ind_nz]
        real_y = real_y[ind_nz]
      else: nt = nt0
      for ni in range(nt):
        layers[s, real_x[ni],real_y[ni]]+=1.0

  # STEP 2  - Convolve Layers and put in pixels
  radius = 1.1
  sig = fwhm / 2.355 / pix
  flattened_pixmap = np.sum(layers,axis=0)
  total_circles_mask = circle_mask(flattened_pixmap, radius * fwhm, pix)
  ind_fit = np.where(total_circles_mask >= 1) # & zeromask != 0)
  nhits = np.shape(ind_fit)[1]
  cfits_maps = np.zeros([nlists,nhits])

  kern = gauss_kern(fwhm, np.floor(fwhm * 10), pix)
  for u in range(nlists):
    layer = layers[u,:,:]
    tmap = smooth_psf(layer, kern)
    tmap[ind_fit] -= np.mean(tmap[ind_fit])
    cfits_maps[u,:] = tmap[ind_fit]

  # STEP 3 - Regress Layers with Map (i.e., stack!)

  cmap[ind_fit] -= np.mean(cmap[ind_fit], dtype=np.float32)

  fit_params = Parameters()

  for iarg in range(nlists):
    fit_params.add('layer'+str(iarg),value= 1e-3*np.random.randn())
  imap = cmap[ind_fit]
  ierr = cnoise[ind_fit]

  cov_ss_1d = minimize(simultaneous_stack_array_oned, fit_params,
    args=(np.ndarray.flatten(cfits_maps),), kws={'data1d':np.ndarray.flatten(imap),'err1d':np.ndarray.flatten(ierr)})

  return cov_ss_1d

def stack_libraries_in_redshift_slices_old(

  map_library,
  subcatalog_library,
  quiet=None):

  #FIRST DO LAYERS CUBE
  n_sources_max=500000l
  nwv = len(map_library.keys())
  lists = subcatalog_library.keys()
  nlists = len(lists)
  stacked_sed=np.zeros([nwv, nlists])
  stacked_sed_err=np.zeros([nwv,nlists])
  stacked_layers = {}

  #PUT DATA INTO CUBE
  nsources = 0 # initialize a counter
  layers_radec = np.zeros([n_sources_max, nlists, 2]) # nsources by nlis/nts by 2 for RA/DEC
  for i in range(nlists):
    subcatalog_key = subcatalog_library.keys()[i]
    if len(subcatalog_library[subcatalog_key][0]) > 0:
      ra  = subcatalog_library[subcatalog_key][0]
      dec = subcatalog_library[subcatalog_key][1]
      nsources_list=len(ra)
      if nsources_list > n_sources_max:
        print str(nsources_list)+' is too many sources in catalog: use N_SOURCES_MAX flag'
        break
      if nsources_list > 0:
        layers_radec[0:nsources_list,i,0]=ra
        layers_radec[0:nsources_list,i,1]=dec
      if nsources_list > nsources:
        nsources=nsources_list
  layers_radec=layers_radec[0:nsources,:,:] # Crop it down to the length of longest list

  cwavelengths = []
  radius = 1.1
  for iwv in range(nwv):
    print 'stacking '+map_library.keys()[iwv]
    #READ MAPS
    cmap = map_library[map_library.keys()[iwv]].map
    cnoise = map_library[map_library.keys()[iwv]].noise
    cwv = map_library[map_library.keys()[iwv]].wavelength
    cwavelengths.append(cwv)
    chd = map_library[map_library.keys()[iwv]].header
    pixsize = map_library[map_library.keys()[iwv]].pixel_size
    kern = map_library[map_library.keys()[iwv]].psf
    fwhm = map_library[map_library.keys()[iwv]].fwhm
    cw = WCS(chd)
    cms = np.shape(cmap)

    # STEP 1  - Make Layers Cube at each wavelength
    layers=np.zeros([nlists,cms[0],cms[1]])

    for s in range(nlists):
      ind_src = np.where(layers_radec[:,s,0] != 0)
      if np.shape(ind_src)[1] > 0:
        ra = layers_radec[ind_src,s,0]
        dec = layers_radec[ind_src,s,1]
        #pdb.set_trace()
        ty,tx = cw.wcs_world2pix(ra, dec, 0)
        # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
        ind_keep = np.where((tx[0] >= 0) & (np.round(tx[0]) < cms[0]) & (ty[0] >= 0) & (np.round(ty[0]) < cms[1]))
        real_x=np.round(tx[0,ind_keep][0]).astype(int)
        real_y=np.round(ty[0,ind_keep][0]).astype(int)
        # CHECK FOR SOURCES THAT FALL ON ZEROS
        ind_nz=np.where(cmap[real_x,real_y] != 0 )
        nt = np.shape(ind_nz)[1]
        #print 'ngals' + str(nt)
        if nt > 0:
          real_x = real_x[ind_nz]
          real_y = real_y[ind_nz]

          for ni in range(nt):
            layers[s, real_x[ni],real_y[ni]]+=1.0

    # STEP 2  - Convolve Layers and put in pixels
    flattened_pixmap = np.sum(layers,axis=0)
    total_circles_mask = circle_mask(flattened_pixmap, radius * fwhm, pixsize)
    ind_fit = np.where(total_circles_mask >= 1) # & zeromask != 0)
    nhits = np.shape(ind_fit)[1]
    cfits_maps = np.zeros([nlists,nhits])
    cfits_flat = np.asarray([])

    #print cms
    #pdb.set_trace()
    for u in range(nlists):
      layer = layers[u,:,:]
      #tmap = pad_and_smooth_psf(layer, kern)
      tmap = smooth_psf(layer, kern)
      tmap[ind_fit] -= np.mean(tmap[ind_fit])
      cfits_flat = np.append(cfits_flat,np.ndarray.flatten(tmap[ind_fit]))
      cfits_maps[u,:] = tmap[ind_fit]

    #print str(cwv)+' cube smoothed'

    cmap[ind_fit] -= np.mean(cmap[ind_fit], dtype=np.float32)
    flat_map = np.ndarray.flatten(cmap[ind_fit])
    flat_noise = np.ndarray.flatten(cnoise[ind_fit])

    fit_params = Parameters()
    for iarg in range(nlists):

      #fit_params.add('layer'+str(iarg),value= 1e-3*np.random.randn())
      #fit_params.add('z_0.5-1.0__m_11.0-13.0_qt'+str(iarg),value= 1e-3*np.random.randn())
      #fit_params.add(lists[iarg],value= 1e-3*np.random.randn())
      arg = clean_args(lists[iarg])
      #arg=arg.replace('.','p')
      #arg=arg.replace('-','_')
      #print arg
      #arg = arg.translate(None, ''.join(['-','.']))
      fit_params.add(arg,value= 1e-3*np.random.randn())
    imap = cmap[ind_fit]
    ierr = cnoise[ind_fit]

    cov_ss_1d = minimize(simultaneous_stack_array_oned, fit_params,
      args=(np.ndarray.flatten(cfits_maps),), kws={'data1d':np.ndarray.flatten(imap),'err1d':np.ndarray.flatten(ierr)})

    stacked_flux = np.array(cov_ss_1d.params.values())
    stacked_sed[iwv,:] = stacked_flux
    stacked_layers[str(cwv)] = cov_ss_1d.params

    #print  map_library.keys()[iwv]+' stack completed'
    #pdb.set_trace()

  ind_sorted = np.argsort(np.asarray(cwavelengths))
  new_stacked_sed = np.array([stacked_sed[i,:] for i in ind_sorted])

  return stacked_layers
  #return new_stacked_sed

def stack_libraries_in_redshift_slices(
  map_library,
  subcatalog_library,
  quiet=None):
  #more tests

  map_names = [i for i in map_library.keys()]
  # All wavelengths in cwavelengths
  cwavelengths = [map_library[i].wavelength for i in map_names]
  # Unique wavelengths in uwavelengths
  uwavelengths = np.sort(np.unique(cwavelengths))
  # nwv the number of unique wavelengths
  nwv = len(uwavelengths)

  lists = subcatalog_library.keys()
  nlists = len(lists)
  stacked_sed=np.zeros([nwv, nlists])
  stacked_sed_err=np.zeros([nwv,nlists])
  stacked_layers = {}

  cwavelengths = []
  radius = 1.1
  for iwv in range(nwv):
    print 'stacking '+map_library.keys()[iwv]
    #READ MAPS
    cmap = map_library[map_library.keys()[iwv]].map
    cnoise = map_library[map_library.keys()[iwv]].noise
    cwv = map_library[map_library.keys()[iwv]].wavelength
    cwavelengths.append(cwv)
    chd = map_library[map_library.keys()[iwv]].header
    pixsize = map_library[map_library.keys()[iwv]].pixel_size
    kern = map_library[map_library.keys()[iwv]].psf
    fwhm = map_library[map_library.keys()[iwv]].fwhm
    cw = WCS(chd)
    cms = np.shape(cmap)

    # STEP 1  - Make Layers Cube at each wavelength
    layers=np.zeros([nlists,cms[0],cms[1]])

    for k in range(nlists):
      s = lists[k]
      if len(subcatalog_library[s][0]) > 0:
        ra = subcatalog_library[s][0]
        dec = subcatalog_library[s][1]
        ty,tx = cw.wcs_world2pix(ra, dec, 0)
        # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
        ind_keep = np.where((np.round(tx) >= 0) & (np.round(tx) < cms[0]) & (np.round(ty) >= 0) & (np.round(ty) < cms[1]))
        real_x=np.round(tx[ind_keep]).astype(int)
        real_y=np.round(ty[ind_keep]).astype(int)
        # CHECK FOR SOURCES THAT FALL ON ZEROS
        ind_nz=np.where(cmap[real_x,real_y] != 0 )
        nt = np.shape(ind_nz)[1]
        #print 'ngals: ' + str(nt)
        if nt > 0:
          real_x = real_x[ind_nz]
          real_y = real_y[ind_nz]

          for ni in range(nt):
            layers[k, real_x[ni],real_y[ni]]+=1.0

    # STEP 2  - Convolve Layers and put in pixels
    flattened_pixmap = np.sum(layers,axis=0)
    total_circles_mask = circle_mask(flattened_pixmap, radius * fwhm, pixsize)
    #pdb.set_trace()
    ind_fit = np.where(total_circles_mask >= 1) # & zeromask != 0)
    nhits = np.shape(ind_fit)[1]
    ###
    #cfits_maps = np.zeros([nlists,nhits])
    cfits_flat = np.asarray([])
    ###

    #print cms
    #pdb.set_trace()
    for u in range(nlists):
      layer = layers[u,:,:]
      #tmap = pad_and_smooth_psf(layer, kern)
      tmap = smooth_psf(layer, kern)
      tmap[ind_fit] -= np.mean(tmap[ind_fit])
      cfits_flat = np.append(cfits_flat,np.ndarray.flatten(tmap[ind_fit]))
      #cfits_maps[u,:] = tmap[ind_fit]

    #print str(cwv)+' cube smoothed'

    cmap[ind_fit] -= np.mean(cmap[ind_fit], dtype=np.float32)
    #flat_map = np.ndarray.flatten(cmap[ind_fit])
    #flat_noise = np.ndarray.flatten(cnoise[ind_fit])
    imap = np.ndarray.flatten(cmap[ind_fit])
    ierr = np.ndarray.flatten(cnoise[ind_fit])

    fit_params = Parameters()
    for iarg in range(nlists):
      #fit_params.add('layer'+str(iarg),value= 1e-3*np.random.randn())
      #fit_params.add('z_0.5-1.0__m_11.0-13.0_qt'+str(iarg),value= 1e-3*np.random.randn())
      #fit_params.add(lists[iarg],value= 1e-3*np.random.randn())
      arg = clean_args(lists[iarg])
      #arg=arg.replace('.','p')
      #arg=arg.replace('-','_')
      #print arg
      #arg = arg.translate(None, ''.join(['-','.']))
      fit_params.add(arg,value= 1e-3*np.random.randn())

    #pdb.set_trace()
    cov_ss_1d = minimize(simultaneous_stack_array_oned, fit_params,
      args=(cfits_flat,), kws={'data1d':imap,'err1d':ierr}, nan_policy = 'propagate')

    stacked_flux = np.array(cov_ss_1d.params.values())
    stacked_sed[iwv,:] = stacked_flux
    stacked_layers[str(cwv)] = cov_ss_1d.params

    #print  map_library.keys()[iwv]+' stack completed'
    #pdb.set_trace()

  #pdb.set_trace()
  ind_sorted = np.argsort(np.asarray(cwavelengths))
  new_stacked_sed = np.array([stacked_sed[i,:] for i in ind_sorted])

  return stacked_layers
  #return new_stacked_sed

def stack_libraries_in_layers(
  map_library,
  subcatalog_library,
  quiet=None):

  map_names = [i for i in map_library.keys()]
  # All wavelengths in cwavelengths
  cwavelengths = [map_library[i].wavelength for i in map_names]
  # Unique wavelengths in uwavelengths
  uwavelengths = np.sort(np.unique(cwavelengths))
  # nwv the number of unique wavelengths
  nwv = len(uwavelengths)

  lists = subcatalog_library.keys()
  nlists = len(lists)
  #stacked_sed=np.zeros([nwv, nlists])
  #stacked_sed_err=np.zeros([nwv,nlists])
  stacked_layers = {}

  cwavelengths = []
  radius = 1.1
  for iwv in range(nwv):
    print 'stacking '+map_library.keys()[iwv]
    #READ MAPS
    cmap = map_library[map_library.keys()[iwv]].map
    cnoise = map_library[map_library.keys()[iwv]].noise
    cwv = map_library[map_library.keys()[iwv]].wavelength
    cname = map_library.keys()[iwv]
    cwavelengths.append(cwv)
    chd = map_library[map_library.keys()[iwv]].header
    pixsize = map_library[map_library.keys()[iwv]].pixel_size
    kern = map_library[map_library.keys()[iwv]].psf
    fwhm = map_library[map_library.keys()[iwv]].fwhm
    cw = WCS(chd)
    cms = np.shape(cmap)

    # STEP 1  - Make Layers Cube at each wavelength
    layers=np.zeros([nlists,cms[0],cms[1]])

    for k in range(nlists):
      s = lists[k]
      if len(subcatalog_library[s][0]) > 0:
        ra = subcatalog_library[s][0]
        dec = subcatalog_library[s][1]
        ty,tx = cw.wcs_world2pix(ra, dec, 0)
        # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
        ind_keep = np.where((np.round(tx) >= 0) & (np.round(tx) < cms[0]) & (np.round(ty) >= 0) & (np.round(ty) < cms[1]))
        real_x=np.round(tx[ind_keep]).astype(int)
        real_y=np.round(ty[ind_keep]).astype(int)
        # CHECK FOR SOURCES THAT FALL ON ZEROS
        ind_nz=np.where(cmap[real_x,real_y] != 0 )
        nt = np.shape(ind_nz)[1]
        #print 'ngals: ' + str(nt)
        if nt > 0:
          real_x = real_x[ind_nz]
          real_y = real_y[ind_nz]

          for ni in range(nt):
            layers[k, real_x[ni],real_y[ni]]+=1.0

    # STEP 2  - Convolve Layers and put in pixels
    flattened_pixmap = np.sum(layers,axis=0)
    total_circles_mask = circle_mask(flattened_pixmap, radius * fwhm, pixsize)
    ind_fit = np.where(total_circles_mask >= 1) # & zeromask != 0)
    del total_circles_mask
    nhits = np.shape(ind_fit)[1]
    ###
    #cfits_maps = np.zeros([nlists,nhits])
    cfits_flat = np.asarray([])
    ###

    #print cms
    #pdb.set_trace()
    for u in range(nlists):
      layer = layers[u,:,:]
      #tmap = pad_and_smooth_psf(layer, kern)
      tmap = smooth_psf(layer, kern)
      tmap[ind_fit] -= np.mean(tmap[ind_fit])
      cfits_flat = np.append(cfits_flat,np.ndarray.flatten(tmap[ind_fit]))
      #cfits_maps[u,:] = tmap[ind_fit]

    #print str(cwv)+' cube smoothed'

    cmap[ind_fit] -= np.mean(cmap[ind_fit], dtype=np.float32)
    #flat_map = np.ndarray.flatten(cmap[ind_fit])
    #flat_noise = np.ndarray.flatten(cnoise[ind_fit])
    imap = np.ndarray.flatten(cmap[ind_fit])
    ierr = np.ndarray.flatten(cnoise[ind_fit])

    fit_params = Parameters()
    for iarg in range(nlists):
      arg = clean_args(lists[iarg])
      fit_params.add(arg,value= 1e-3*np.random.randn())


    if len(ierr)==0: pdb.set_trace()
    cov_ss_1d = minimize(simultaneous_stack_array_oned, fit_params,
      args=(cfits_flat,), kws={'data1d':imap,'err1d':ierr}, nan_policy = 'propagate')
    del cfits_flat, imap, ierr

    #Dictionary keys decided here.  Was originally wavelengths.  Changing it back to map_names
    #stacked_layers[str(cwv)] = cov_ss_1d.params
    stacked_layers[cname] = cov_ss_1d.params

    #print  map_library.keys()[iwv]+' stack completed'
    #pdb.set_trace()

  gc.collect
  return stacked_layers

def stack_multiple_fields_in_redshift_slices_old(
  map_library,
  subcatalog_library,
  quiet=None):

  n_sources_max=500000l
  map_names = [i for i in map_library.keys()]
  # All wavelengths in cwavelengths
  cwavelengths = [map_library[i].wavelength for i in map_names]
  # Unique wavelengths in uwavelengths
  uwavelengths = np.sort(np.unique(cwavelengths))
  # nwv the number of unique wavelengths
  nwv = len(uwavelengths)
  #pdb.set_trace()

  lists = subcatalog_library.keys()
  nlists = len(lists)
  stacked_sed=np.zeros([nwv, nlists])
  stacked_sed_err=np.zeros([nwv,nlists])
  stacked_layers = {}

  #PUT DATA INTO CUBE
  nsources = 0 # initialize a counter
  #THIS IS BAD FORM, SHOULD BE TURNED INTO DICTIONARIES
  layers_radec = np.zeros([n_sources_max, nlists, 2]) # nsources by nlis/nts by 2 for RA/DEC
  #THIS SEEMS LIKE ITS TROUBLE....
  for i in range(nlists):
    subcatalog_key = subcatalog_library.keys()[i]
    if len(subcatalog_library[subcatalog_key][0]) > 0:
      ra  = subcatalog_library[subcatalog_key][0]
      dec = subcatalog_library[subcatalog_key][1]
      nsources_list=len(ra)
      if nsources_list > n_sources_max:
        print str(nsources_list)+' is too many sources in catalog: use N_SOURCES_MAX flag'
        break
      if nsources_list > 0:
        layers_radec[0:nsources_list,i,0]=ra
        layers_radec[0:nsources_list,i,1]=dec
      if nsources_list > nsources:
        nsources=nsources_list
  layers_radec=layers_radec[0:nsources,:,:] # Crop it down to the length of longest list

  radius = 1.1
  #pdb.set_trace()
  for jwv in range(nwv):
    argwv = np.where(cwavelengths == uwavelengths[jwv])[0]
    ninstances = cwavelengths.count(uwavelengths[jwv])
    for iwv in argwv:
      print map_names[iwv]
      print 'stacking '+str(ninstances)+' maps at ' + str(cwavelengths[iwv])
      #pdb.set_trace()
      #NOW NEED TO FIGURE HOW HOW TO LOOP BELOW
      cfits_flat = np.asarray([])
      imap = np.asarray([])
      ierr = np.asarray([])
      #READ MAPS
      cmap = map_library[map_names[iwv]].map
      cnoise = map_library[map_names[iwv]].noise
      #cwv = map_library[map_names[iwv]].wavelength
      #cwavelengths.append(cwv)
      chd = map_library[map_names[iwv]].header
      pixsize = map_library[map_names[iwv]].pixel_size
      kern = map_library[map_names[iwv]].psf
      fwhm = map_library[map_names[iwv]].fwhm
      cw = WCS(chd)
      cms = np.shape(cmap)

      # STEP 1  - Make Layers Cube at each wavelength
      layers=np.zeros([nlists,cms[0],cms[1]])

      for s in range(nlists):
        ind_src = np.where(layers_radec[:,s,0] != 0)
        if np.shape(ind_src)[1] > 0:
          ra = layers_radec[ind_src,s,0]
          dec = layers_radec[ind_src,s,1]
          ty,tx = cw.wcs_world2pix(ra, dec, 0)
          # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
          ind_keep = np.where((tx[0] >= 0) & (np.round(tx[0]) < cms[0]) & (ty[0] >= 0) & (np.round(ty[0]) < cms[1]))
          real_x=np.round(tx[0,ind_keep][0]).astype(int)
          real_y=np.round(ty[0,ind_keep][0]).astype(int)
          # CHECK FOR SOURCES THAT FALL ON ZEROS
          ind_nz=np.where(cmap[real_x,real_y] != 0 )
          nt = np.shape(ind_nz)[1]
          #print 'ngals' + str(nt)
          if nt > 0:
            real_x = real_x[ind_nz]
            real_y = real_y[ind_nz]

            for ni in range(nt):
              layers[s, real_x[ni],real_y[ni]]+=1.0

      # STEP 2  - Convolve Layers and put in pixels
      flattened_pixmap = np.sum(layers,axis=0)
      total_circles_mask = circle_mask(flattened_pixmap, radius * fwhm, pixsize)
      ind_fit = np.where(total_circles_mask >= 1) # & zeromask != 0)
      nhits = np.shape(ind_fit)[1]
      #cfits_maps = np.zeros([nlists,nhits])

      #print cms
      #pdb.set_trace()
      for u in range(nlists):
        layer = layers[u,:,:]
        tmap = pad_and_smooth_psf(layer, kern)
        tmap[ind_fit] -= np.mean(tmap[ind_fit])
        cfits_flat = np.append(cfits_flat,np.ndarray.flatten(tmap[ind_fit]))
        #cfits_maps[u,:] = tmap[ind_fit]

      #print str(cwv)+' cube smoothed'

      cmap[ind_fit] -= np.mean(cmap[ind_fit], dtype=np.float32)
      flat_map = np.ndarray.flatten(cmap[ind_fit])
      flat_noise = np.ndarray.flatten(cnoise[ind_fit])

      fit_params = Parameters()
      for iarg in range(nlists):
        arg = clean_args(lists[iarg])
        #arg=arg.replace('.','p')
        #arg=arg.replace('-','_')
        fit_params.add(arg,value= 1e-3*np.random.randn())
      imap = np.append(imap,np.ndarray.flatten(cmap[ind_fit]))
      ierr = np.append(ierr,np.ndarray.flatten(cnoise[ind_fit]))

    #END MULTIMAP LOOP AND STACK HERE

    cov_ss_1d = minimize(simultaneous_stack_array_oned, fit_params,
      args=(cfits_flat,), kws={'data1d':imap,'err1d':ierr})
      #args=(np.ndarray.flatten(cfits_maps),), kws={'data1d':np.ndarray.flatten(imap),'err1d':np.ndarray.flatten(ierr)})

    stacked_flux = np.array(cov_ss_1d.params.values())
    stacked_sed[jwv,:] = stacked_flux
    stacked_layers[str(uwavelengths[jwv])] = cov_ss_1d.params

    #print  map_library.keys()[iwv]+' stack completed'
    #pdb.set_trace()

  ind_sorted = np.argsort(np.asarray(uwavelengths))
  new_stacked_sed = np.array([stacked_sed[i,:] for i in ind_sorted])

  return stacked_layers
  #return new_stacked_sed

def stack_multiple_fields_in_redshift_slices(
  map_library,
  subcatalog_library,
  quiet=None):

  map_names = [i for i in map_library.keys()]
  # All wavelengths in cwavelengths
  cwavelengths = [map_library[i].wavelength for i in map_names]
  # Unique wavelengths in uwavelengths
  uwavelengths = np.sort(np.unique(cwavelengths))
  # nwv the number of unique wavelengths
  nwv = len(uwavelengths)

  lists = subcatalog_library.keys()
  nlists = len(lists)
  stacked_sed=np.zeros([nwv, nlists])
  stacked_sed_err=np.zeros([nwv,nlists])
  stacked_layers = {}

  radius = 1.1
  for jwv in range(nwv):
    argwv = np.where(cwavelengths == uwavelengths[jwv])[0]
    ninstances = cwavelengths.count(uwavelengths[jwv])
    cfits_flat = np.asarray([])
    imap = np.asarray([])
    ierr = np.asarray([])
    lmaps = []
    print 'stacking '+str(ninstances)+' maps at ' + str(uwavelengths[jwv])
    for iwv in argwv:
      print map_names[iwv]
      #pdb.set_trace()
      #READ MAPS
      cmap = map_library[map_names[iwv]].map
      cnoise = map_library[map_names[iwv]].noise
      chd = map_library[map_names[iwv]].header
      pixsize = map_library[map_names[iwv]].pixel_size
      kern = map_library[map_names[iwv]].psf
      fwhm = map_library[map_names[iwv]].fwhm
      cw = WCS(chd)
      cms = np.shape(cmap)

      # STEP 1  - Make Layers Cube for each map
      layers=np.zeros([nlists,cms[0],cms[1]])

      for k in range(nlists):
        s = lists[k]
        if len(subcatalog_library[s][0]) > 0:
          ra = subcatalog_library[s][0]
          dec = subcatalog_library[s][1]
          #
          ty,tx = cw.wcs_world2pix(ra, dec, 0)
          # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
          ind_keep = np.where((np.round(tx) >= 0) & (np.round(tx) < cms[0]) & (np.round(ty) >= 0) & (np.round(ty) < cms[1]))
          real_x=np.round(tx[ind_keep]).astype(int)
          real_y=np.round(ty[ind_keep]).astype(int)
          # CHECK FOR SOURCES THAT FALL ON ZEROS
          ind_nz=np.where(cmap[real_x,real_y] != 0 )
          nt = np.shape(ind_nz)[1]
          #print 'ngals: ' + str(nt)
          if nt > 0:
            real_x = real_x[ind_nz]
            real_y = real_y[ind_nz]
            for ni in range(nt):
              layers[k, real_x[ni],real_y[ni]]+=1.0

      # STEP 2  - Convolve Layers and put in pixels
      flattened_pixmap = np.sum(layers,axis=0)
      total_circles_mask = circle_mask(flattened_pixmap, radius * fwhm, pixsize)
      ind_fit = np.where(total_circles_mask >= 1) # & zeromask != 0)
      nhits = np.shape(ind_fit)[1]

      #print 'there'
      for u in range(nlists):
        #print lists[u]
        layer = layers[u,:,:]
        #tmap = pad_and_smooth_psf(layer, kern)
        tmap = smooth_psf(layer, kern)
        tmap[ind_fit] -= np.mean(tmap[ind_fit], dtype=np.float32)
        cfits_flat = np.append(cfits_flat,np.ndarray.flatten(tmap[ind_fit]))

      #print str(cwv)+' cube smoothed'

      cmap[ind_fit] -= np.mean(cmap[ind_fit], dtype=np.float32)
      flat_map = np.ndarray.flatten(cmap[ind_fit])
      flat_noise = np.ndarray.flatten(cnoise[ind_fit])

      imap = np.append(imap,flat_map) # np.ndarray.flatten(cmap[ind_fit]))
      ierr = np.append(ierr,flat_noise) # np.ndarray.flatten(cnoise[ind_fit]))
      #print 'imap length: ' + str(len(imap))
      #print 'ierr length: ' + str(len(ierr))
      lmaps = lmaps + [len(flat_map)]

    #END MULTIMAP LOOP AND STACK HERE
    fit_params = Parameters()
    arg_order = []
    for arg in lists:
      #arg=arg.replace('.','p')
      #arg=arg.replace('-','_')
      arg_order.append(clean_args(arg))
      fit_params.add(arg,value= 1e-3*np.random.randn())

    #pdb.set_trace()
    print 'cfits_flat length: ' + str(len(cfits_flat))
    cov_ss_1d = minimize(simultaneous_stack_multimap, fit_params,
      args=(cfits_flat,), kws={'data1d':imap,'err1d':ierr,'nmaps':ninstances, 'lmaps':lmaps})

    stacked_flux = np.array(cov_ss_1d.params.values())
    stacked_sed[jwv,:] = stacked_flux
    stacked_layers[str(uwavelengths[jwv])] = cov_ss_1d.params

    #print  map_library.keys()[iwv]+' stack completed'

  ind_sorted = np.argsort(np.asarray(uwavelengths))
  new_stacked_sed = np.array([stacked_sed[i,:] for i in ind_sorted])

  return stacked_layers
  #return new_stacked_sed
