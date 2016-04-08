import pdb
import numpy as np
from astropy.wcs import WCS
from utils import circle_mask
from utils import dist_idl
from utils import gauss_kern
from utils import pad_and_smooth_psf
from utils import shift_twod
from utils import smooth_psf
from lmfit import Parameters, minimize, fit_report

def simultaneous_stack_array_oned(p, layers_1d, data1d, err1d = None):
  ''' Function to Minimize written specifically for lmfit '''

  v = p.valuesdict()

  len_model = len(data1d)
  nlayers = len(layers_1d)/len_model

  model = np.zeros(len_model)

  for i in range(nlayers):
    #model[:] += layers_1d[i*len_model:(i+1)*len_model] * v['layer'+str(i)] 
    model[:] += layers_1d[i*len_model:(i+1)*len_model] * v[v.keys()[i]] 

  if err1d is None:
    return (data1d - model)
  return (data1d - model)/err1d

def simultaneous_stack_array(p, layers_2d, data, err = None):
  ''' Function to Minimize written specifically for lmfit '''

  v = p.valuesdict()

  csize = np.shape(layers_2d)

  model = np.zeros(csize[1])

  for i in range(csize[0]):
    #model += layers_2d[i,:] * v['layer'+str(i)] 
    model += layers_2d[i,:] * v[v.keys()[i]]  

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

def stack_libraries_in_redshift_slices(
  map_library, 
  subcatalog_library,
  quiet=None):
  
  #FIRST DO LAYERS CUBE
  n_sources_max=50000l
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
        print 'too many sources in catalog: use N_SOURCES_MAX flag'
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
      tmap = pad_and_smooth_psf(layer, kern)
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
      arg = lists[iarg]
      arg=arg.replace('.','p')
      arg=arg.replace('-','_')
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

def stack_multiple_fields_in_redshift_slices(
  map_library, 
  subcatalog_library,
  quiet=None):
  
  n_sources_max=50000l
  ckeys = [i for i in map_library.keys()]
  cwavelengths = [map_library[i].wavelength for i in ckeys] 
  #indwv=np.argsort(cwavelengths)
  uwavelengths = np.sort(np.unique(cwavelengths))
  nwv = len(uwavelengths)
  #pdb.set_trace()

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
        print 'too many sources in catalog: use N_SOURCES_MAX flag'
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
      print ckeys[iwv]
      print 'stacking '+str(ninstances)+' maps at ' + str(cwavelengths[iwv]) 
      #pdb.set_trace()
      #NOW NEED TO FIGURE HOW HOW TO LOOP BELOW
      cfits_flat = np.asarray([])
      imap = np.asarray([])
      ierr = np.asarray([])
      #READ MAPS
      cmap = map_library[ckeys[iwv]].map
      cnoise = map_library[ckeys[iwv]].noise
      #cwv = map_library[ckeys[iwv]].wavelength
      #cwavelengths.append(cwv)
      chd = map_library[ckeys[iwv]].header
      pixsize = map_library[ckeys[iwv]].pixel_size
      kern = map_library[ckeys[iwv]].psf
      fwhm = map_library[ckeys[iwv]].fwhm
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
        arg = lists[iarg]
        arg=arg.replace('.','p')
        arg=arg.replace('-','_')
        fit_params.add(arg,value= 1e-3*np.random.randn())
      imap = np.append(imap,np.ndarray.flatten(cmap[ind_fit]))
      ierr = np.append(ierr,np.ndarray.flatten(cnoise[ind_fit]))

    #END MULTIMAP LOOP AND STACK HERE

    cov_ss_1d = minimize(simultaneous_stack_array_oned, fit_params, 
      args=(cfits_flat,), kws={'data1d':imap,'err1d':ierr})
      #args=(np.ndarray.flatten(cfits_maps),), kws={'data1d':np.ndarray.flatten(imap),'err1d':np.ndarray.flatten(ierr)})

    stacked_flux = np.array(cov_ss_1d.params.values())
    stacked_sed[jwv,:] = stacked_flux 
    stacked_layers[uwavelengths[jwv]] = cov_ss_1d.params

    #print  map_library.keys()[iwv]+' stack completed'
    #pdb.set_trace()

  ind_sorted = np.argsort(np.asarray(uwavelengths))
  new_stacked_sed = np.array([stacked_sed[i,:] for i in ind_sorted])

  return stacked_layers
  #return new_stacked_sed





