import pdb
import numpy as np
import os
from astropy.wcs import WCS
from invert_sed import single_simple_flux_from_greybody
import readcol
from utils import circle_mask
from utils import dist_idl
from utils import gauss_kern
from utils import pad_and_smooth_psf
from utils import shift_twod
from utils import smooth_psf
from utils import zero_pad
from lmfit import Parameters, minimize, fit_report

def simultaneous_stack_sed_oned(p, layers_1d, data1d, wavelengths, LenLayers, zed, err1d = None):
  ''' Function to Minimize written specifically for lmfit '''

  v = p.valuesdict()

  nwv = len(wavelengths)
  LenModel = len(data1d) 
  Nlayers = len(layers_1d)/LenModel

  model = np.zeros(LenModel)

  LL = [0]
  SuperChunk = [0]
  LL.extend(LenLayers)
  cuLL = np.cumsum(LL)
  SuperChunk.extend(LenLayers * Nlayers)
  cuSuperChunk = np.cumsum(SuperChunk)
  pdb.set_trace()
  for i in range(Nlayers):
    Temp = np.asarray(v['T'+str(i)])
    Lir = np.asarray(v['L'+str(i)])

    fluxes = single_simple_flux_from_greybody(np.asarray(wavelengths), Trf = Temp, Lrf = Lir, b=2, zin=zed)
    for iwv in range(nwv):
      model[cuLL[iwv]:cuLL[iwv+1]] +=  fluxes[0][iwv] * layers_1d[ cuSuperChunk[iwv] + i * LenLayers[iwv]: cuSuperChunk[iwv] + (i+1) * LenLayers[iwv] ] 

  if err1d is None:
    return (data1d - model)
  return (data1d - model)/err1d

def stack_in_redshift_slices(
  cmaps, 
  hd, 
  layers_radec, 
  wavelengths,
  fwhm=None, 
  psf_names=None,
  cnoise=None, 
  mask=None, 
  beam_area=None, 
  err_ss=None, 
  zed=0.01,
  quiet=None):
  
  w = WCS(hd)
  #FIND SIZES OF MAP AND LISTS
  cms = np.shape(cmaps) # should be a cube
  nwv = cms[0] 
  #zeromask = np.zeros(cms)

  size_cube = np.shape(layers_radec)
  nsrcmax = size_cube[0]
  nlists = int(size_cube[1])
  
  ind_map_zero = np.where(np.isnan(cmaps))
  nzero = np.shape(ind_map_zero)[1]

  if np.sum(cnoise) == 0: cnoise=cmaps*0.0 + 1.0

  pix=hd["CD2_2"]*3600.
  if pix == 0: pix=hd["CDELT2"]*3600.

  for iwv in range(nwv):
    #[STEP 0] - Calibrate maps
    if beam_area != None:
      cmaps[iwv,:,:]=cmaps[iwv,:,:]*beam_area[iwv]*1e6
      cnoise[iwv,:,:]=noise[iwv,:,:]*beam_area[iwv]*1e6

  # STEP 1  - Make Layers Cube
  layers=np.zeros([nlists,cms[1],cms[2]])

  for s in range(nlists):
    ind_src = np.where(layers_radec[:,s,0] != 0)
    if np.shape(ind_src)[1] > 0:
      ra = layers_radec[ind_src,s,0]
      dec = layers_radec[ind_src,s,1]
      ty,tx = w.wcs_world2pix(ra, dec, 0) 
      # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
      ind_keep = np.where((tx[0] >= 0) & (np.round(tx[0]) < cms[1]) & (ty[0] >= 0) & (np.round(ty[0]) < cms[2]))
      nt0 = np.shape(ind_keep)[1]
      real_x=np.round(tx[0,ind_keep][0]).astype(int)
      real_y=np.round(ty[0,ind_keep][0]).astype(int)
      # CHECK FOR SOURCES THAT FALL ON ZEROS MAP
      # THIS NEEDS COMPLETE OVERHAUL, PARTICULARLY WHEN INCLUDING DIFFERENT AREA MAPS!!
      if nzero > 0:
        tally = np.zeros(nt0)
        for d in range(nt0):
          if cmaps[0,real_x[d],real_y[d]] != 0: 
            tally[d]=1.
        ind_nz=np.where(tally == 1)
        nt = np.shape(ind_nz)[1]
        real_x = real_x[ind_nz]
        real_y = real_y[ind_nz]
      else: nt = nt0
      for ni in range(nt):
        layers[s, real_x[ni],real_y[ni]]+=1.0

  # STEP 2  - Convolve Layers and put in pixels
  #all_map_layers = np.zeros(np.append(nwv,np.shape(layers)))

  cfits_flat = np.asarray([])
  cfits_flat2= np.asarray([])
  flat_maps= np.asarray([])
  flat_noise= np.asarray([])
  LenLayers= np.zeros([nwv])

  radius = 1.1
  for iwv in range(nwv):
    sig = fwhm[iwv] / 2.355 / pix 
    flattened_pixmap = np.sum(layers,axis=0)
    total_circles_mask = circle_mask(flattened_pixmap, radius * fwhm[iwv], pix)
    ind_fit = np.where(total_circles_mask >= 1) # & zeromask != 0)
    nhits = np.shape(ind_fit)[1]
    LenLayers[iwv] = nhits

    kern = gauss_kern(fwhm[iwv], np.floor(fwhm[iwv] * 10), pix)
    for u in range(nlists):
      layer = layers[u,:,:]  
      tmap = smooth_psf(layer, kern)
      tmap[ind_fit] -= np.mean(tmap[ind_fit])
      cfits_flat = np.append(cfits_flat,np.ndarray.flatten(tmap[ind_fit]))

    lmap = cmaps[iwv]
    lnoise = cnoise[iwv]
    lmap[ind_fit] -= np.mean(lmap[ind_fit], dtype=np.float32)
    flat_maps = np.append(flat_maps,np.ndarray.flatten(lmap[ind_fit]))
    flat_noise = np.append(flat_noise,np.ndarray.flatten(lnoise[ind_fit]))

  # STEP 3 - Regress Layers with Map (i.e., stack!)

  fit_params = Parameters()

  fit_params.add('b',value= 2.0,vary=False)
  for iarg in range(nlists): 
    fit_params.add('T'+str(iarg),value= 25.,vary=True,min=8.,max=80.)
    fit_params.add('L'+str(iarg),value= 1e12,min=0.,max=1e14)

  pdb.set_trace()
  cov_ss_1d = minimize(simultaneous_stack_sed_oned, fit_params, 
    args=(cfits_flat,), kws={'data1d':flat_maps,'err1d':flat_noise,'wavelengths':wavelengths,'LenLayers':LenLayers,'zed':zed})
    
  return cov_ss_1d

def stack_objects_in_redshift_slices(
  map_library, 
  #catalog_library, #would have as keys names of variables...
  subcatalog_names,
  zed=0.01,
  quiet=None):

  #FIRST DO LAYERS CUBE
  n_sources_max=50000l
  nmap = len(map_library.keys())

  #PUT DATA INTO CUBE
  nlists = len(subcatalog_names)
  nsources = 0 # initialize a counter
  layers_radec = np.zeros([n_sources_max, nlists, 2]) # nsources by nlis/nts by 2 for RA/DEC
  for i in range(nlists): 
    list_name = subcatalog_names[i]
    if os.path.getsize(list_name) > 0: 
      ra, dec = readcol.readcol(list_name,fsep=',',twod=False)
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
  #layers_radec=layers_radec[0:nsources-1,:,:] # Crop it down to the length of longest list
  stacked_sed=np.zeros([nmap, nlists])

  cmaps = [] 
  cnoise = [] 
  cwavelengths = []
  cw = []
  cpix = []
  cms = []
  ckern = []
  cfwhm = []
  for wv in range(nmap): 
    print map_library.keys()[wv]
    #READ MAPS
    tmaps = map_library[map_library.keys()[wv]].map
    tnoise = map_library[map_library.keys()[wv]].noise
    #if beam_area != None:
    #  tmaps*=tmaps*beam_area[iwv]*1e6
    #  tnoise*=toise*beam_area[iwv]*1e6
    twv = map_library[map_library.keys()[wv]].wavelength
    hd = map_library[map_library.keys()[wv]].header
    pixsize = map_library[map_library.keys()[wv]].pixel_size
    kern = map_library[map_library.keys()[wv]].psf
    fwhm = map_library[map_library.keys()[wv]].fwhm
    #color_correction = map_library[map_library.keys()[wv]].color_correction
    #tmaps *= color_correction[wv]
    #tnoise *= color_correction[wv]
    cmaps.append(tmaps)
    cnoise.append(tnoise)
    cwavelengths.append(twv)
    cw.append(WCS(hd))
    cpix.append(pixsize)
    cms.append(np.shape(tmaps))
    ckern.append(kern)
    cfwhm.append(fwhm)
  
  #FIND SIZES OF MAP AND LISTS
  nwv = len(cwavelengths)  

  cfits_flat = np.asarray([])
  cfits_flat2= np.asarray([])
  flat_maps= np.asarray([])
  flat_noise= np.asarray([])
  LenLayers= np.zeros([nwv])

  radius = 1.1
  for iwv in range(nwv):
    # STEP 1  - Make Layers Cube at each wavelength
    layers=np.zeros([nlists,cms[iwv][0],cms[iwv][1]])

    x0 = 1e6
    x1 = 0
    y0 = 1e6
    y1 = 0
    for s in range(nlists):
      ind_src = np.where(layers_radec[:,s,0] != 0)
      if np.shape(ind_src)[1] > 0:
        ra = layers_radec[ind_src,s,0]
        dec = layers_radec[ind_src,s,1]
        ty,tx = cw[iwv].wcs_world2pix(ra, dec, 0) 
        # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
        ind_keep = np.where((tx[0] >= 0) & (np.round(tx[0]) < cms[iwv][0]) & (ty[0] >= 0) & (np.round(ty[0]) < cms[iwv][1]))
        real_x=np.round(tx[0,ind_keep][0]).astype(int)
        real_y=np.round(ty[0,ind_keep][0]).astype(int)
        # CHECK FOR SOURCES THAT FALL ON ZEROS 
        ind_nz=np.where(cmaps[iwv][real_x,real_y] != 0 )
        nt = np.shape(ind_nz)[1]
        #print 'ngals' + str(nt)
        if nt > 0:
          real_x = real_x[ind_nz]
          real_y = real_y[ind_nz]
          if np.min(real_x) < x0: x0 = np.min(real_x) 
          if np.max(real_x) > x1: x1 = np.max(real_x) 
          if np.min(real_y) < y0: y0 = np.min(real_y) 
          if np.max(real_y) > y1: y1 = np.max(real_y) 
          for ni in range(nt):
            layers[s, real_x[ni],real_y[ni]]+=1.0
    # STEP 1b  - Crop maps and layers before convolving
    #print x0, x1, y0, y1
    #bpad = np.ceil(radius * cfwhm[iwv] /pixsize)
    #layers = layers[:,x0-bpad:x1+bpad,y0-bpad:y1+bpad]
    #layers2 = layers[:,x0-bpad:x1+bpad,y0-bpad:y1+bpad]

    # STEP 2  - Convolve Layers and put in pixels
    flattened_pixmap = np.sum(layers,axis=0)
    total_circles_mask = circle_mask(flattened_pixmap, radius * cfwhm[iwv], cpix[iwv])
    ind_fit = np.where(total_circles_mask >= 1) # & zeromask != 0)
    nhits = np.shape(ind_fit)[1]
    LenLayers[iwv] = nhits

    #flattened_pixmap2 = np.sum(layers2,axis=0)
    #total_circles_mask2= circle_mask(flattened_pixmap2, radius * cfwhm[iwv], cpix[iwv])
    #ind_fit2 = np.where(total_circles_mask2 >= 1) # & zeromask != 0)
    #nhits2 = np.shape(ind_fit2)[1]
    #LenLayers2[iwv] = nhits2

    for u in range(nlists):
      #print u
      layer = layers[u,:,:]  
      #layer2 = layers2[u,:,:]  
      #tmap = smooth_psf(layer, ckern[iwv])
      #tmap2= smooth_psf(layer2, ckern[iwv])
      tmap = pad_and_smooth_psf(layer, ckern[iwv])
      #tmap2= pad_and_smooth_psf(layer2, ckern[iwv])

      tmap[ind_fit] -= np.mean(tmap[ind_fit])
      cfits_flat = np.append(cfits_flat,np.ndarray.flatten(tmap[ind_fit]))

      #tmap2[ind_fit2] -= np.mean(tmap[ind_fit2])
      #cfits_flat2 = np.append(cfits_flat2,np.ndarray.flatten(tmap2[ind_fit2]))
    print str(cwavelengths[iwv])+' cube smoothed'

    #pdb.set_trace()
    lmap = cmaps[iwv]#[x0-bpad:x1+bpad,y0-bpad:y1+bpad]
    lnoise = cnoise[iwv]#[x0-bpad:x1+bpad,y0-bpad:y1+bpad]
    #lmap2 = cmaps[iwv][x0-bpad:x1+bpad,y0-bpad:y1+bpad]
    #lnoise2 = cnoise[iwv][x0-bpad:x1+bpad,y0-bpad:y1+bpad]
    lmap[ind_fit] -= np.mean(lmap[ind_fit], dtype=np.float32)
    #lmap2[ind_fit2] -= np.mean(lmap2[ind_fit2], dtype=np.float32)
    flat_maps = np.append(flat_maps,np.ndarray.flatten(lmap[ind_fit]))
    flat_noise = np.append(flat_noise,np.ndarray.flatten(lnoise[ind_fit]))

    #pdb.set_trace()

  # STEP 3 - Regress Layers with Map (i.e., stack!)

  fit_params = Parameters()
  fit_params.add('b',value= 2.0,vary=False)
  for iarg in range(nlists): 
    fit_params.add('T'+str(iarg),value= 25.,vary=True,min=8.,max=80.)
    fit_params.add('L'+str(iarg),value= 1e12,min=0.,max=1e14)

  #print np.min(cfits_flat)
  #print np.max(cfits_flat)
  #print np.size(cfits_flat)
  #pdb.set_trace()
  cov_ss_1d = minimize(simultaneous_stack_sed_oned, fit_params, 
    args=(cfits_flat,), kws={'data1d':flat_maps,'err1d':flat_noise,'wavelengths':cwavelengths,'LenLayers':LenLayers,'zed':zed})

  #return cov_ss_1d
  v = cov_ss_1d.params.valuesdict()

  beta = np.asarray(v['b'])
  for ised in range(nlists):
    Temp = np.asarray(v['T'+str(ised)])
    Lir = np.asarray(v['L'+str(ised)])
    #pdb.set_trace()
    stacked_sed[:,ised]=single_simple_flux_from_greybody(np.sort(np.asarray(cwavelengths)), Trf = Temp, Lrf = Lir, b=beta, zin=zed)
  #pdb.set_trace()
  return [stacked_sed, v]


def stack_libraries_in_redshift_slices(
  map_library, 
  subcatalog_library,
  zed=0.01,
  quiet=None):

  #FIRST DO LAYERS CUBE
  n_sources_max=50000l
  nwv = len(map_library.keys())
  lists = subcatalog_library.keys()
  nlists = len(lists)

  #PUT DATA INTO CUBE
  ##########REPLACE WITH LIBRARY
  nsources = 0 # initialize a counter
  layers_radec = np.zeros([n_sources_max, nlists, 2]) # nsources by nlis/nts by 2 for RA/DEC
  for i in range(nlists): 
    #subcatalog_key = subcatalog_library.keys()[i]
    subcatalog_key = lists[i]
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
  #layers_radec=layers_radec[0:nsources-1,:,:] # Crop it down to the length of longest list
  layers_radec=layers_radec[0:nsources,:,:] # Crop it down to the length of longest list
  stacked_sed=np.zeros([nwv, nlists])
  ######################
  cmaps = [] 
  cnoise = [] 
  cwavelengths = []
  cw = []
  cpix = []
  cms = []
  ckern = []
  cfwhm = []
  for wv in range(nwv): 
    print map_library.keys()[wv]
    #READ MAPS
    tmaps = map_library[map_library.keys()[wv]].map
    tnoise = map_library[map_library.keys()[wv]].noise
    twv = map_library[map_library.keys()[wv]].wavelength
    hd = map_library[map_library.keys()[wv]].header
    pixsize = map_library[map_library.keys()[wv]].pixel_size
    kern = map_library[map_library.keys()[wv]].psf
    fwhm = map_library[map_library.keys()[wv]].fwhm
    cmaps.append(tmaps)
    cnoise.append(tnoise)
    cwavelengths.append(twv)
    cw.append(WCS(hd))
    cpix.append(pixsize)
    cms.append(np.shape(tmaps))
    ckern.append(kern)
    cfwhm.append(fwhm)
  
  #FIND SIZES OF MAP AND LISTS
  #nwv = len(cwavelengths)  

  cfits_flat = np.asarray([])
  flat_maps= np.asarray([])
  flat_noise= np.asarray([])
  LenLayers= np.zeros([nwv])

  radius = 1.1
  for iwv in range(nwv):
    # STEP 1  - Make Layers Cube at each wavelength
    layers=np.zeros([nlists,cms[iwv][0],cms[iwv][1]])

    x0 = 1e6
    x1 = 0
    y0 = 1e6
    y1 = 0
    for s in range(nlists):
      ind_src = np.where(layers_radec[:,s,0] != 0)
      if np.shape(ind_src)[1] > 0:
        ra = layers_radec[ind_src,s,0]
        dec = layers_radec[ind_src,s,1]
        ty,tx = cw[iwv].wcs_world2pix(ra, dec, 0) 
        # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
        ind_keep = np.where((tx[0] >= 0) & (np.round(tx[0]) < cms[iwv][0]) & (ty[0] >= 0) & (np.round(ty[0]) < cms[iwv][1]))
        real_x=np.round(tx[0,ind_keep][0]).astype(int)
        real_y=np.round(ty[0,ind_keep][0]).astype(int)
        # CHECK FOR SOURCES THAT FALL ON ZEROS 
        ind_nz=np.where(cmaps[iwv][real_x,real_y] != 0 )
        nt = np.shape(ind_nz)[1]
        #print 'ngals' + str(nt)
        if nt > 0:
          real_x = real_x[ind_nz]
          real_y = real_y[ind_nz]
          if np.min(real_x) < x0: x0 = np.min(real_x) 
          if np.max(real_x) > x1: x1 = np.max(real_x) 
          if np.min(real_y) < y0: y0 = np.min(real_y) 
          if np.max(real_y) > y1: y1 = np.max(real_y) 
          for ni in range(nt):
            layers[s, real_x[ni],real_y[ni]]+=1.0
    # STEP 1b  - Crop maps and layers before convolving
    #print x0, x1, y0, y1
    #bpad = np.ceil(radius * cfwhm[iwv] /pixsize)
    #layers = layers[:,x0-bpad:x1+bpad,y0-bpad:y1+bpad]
    #layers2 = layers[:,x0-bpad:x1+bpad,y0-bpad:y1+bpad]

    # STEP 2  - Convolve Layers and put in pixels
    flattened_pixmap = np.sum(layers,axis=0)
    total_circles_mask = circle_mask(flattened_pixmap, radius * cfwhm[iwv], cpix[iwv])
    ind_fit = np.where(total_circles_mask >= 1) # & zeromask != 0)
    nhits = np.shape(ind_fit)[1]
    LenLayers[iwv] = nhits

    for u in range(nlists):
      layer = layers[u,:,:]  
      tmap = pad_and_smooth_psf(layer, ckern[iwv])

      tmap[ind_fit] -= np.mean(tmap[ind_fit])
      cfits_flat = np.append(cfits_flat,np.ndarray.flatten(tmap[ind_fit]))

    print str(cwavelengths[iwv])+' cube smoothed'

    #pdb.set_trace()
    lmap = cmaps[iwv]#[x0-bpad:x1+bpad,y0-bpad:y1+bpad]
    lnoise = cnoise[iwv]#[x0-bpad:x1+bpad,y0-bpad:y1+bpad]
    lmap[ind_fit] -= np.mean(lmap[ind_fit], dtype=np.float32)
    flat_maps = np.append(flat_maps,np.ndarray.flatten(lmap[ind_fit]))
    flat_noise = np.append(flat_noise,np.ndarray.flatten(lnoise[ind_fit]))
    #pdb.set_trace()

  # STEP 3 - Regress Layers with Map (i.e., stack!)

  tguess = 27.0*((1.+zed)/(1.+1.))**(0.4)
  fit_params = Parameters()
  fit_params.add('b',value= 2.0,vary=False)
  for iarg in range(nlists): 
    arg = lists[iarg]
    arg=arg.replace('.','p')
    arg=arg.replace('-','_')
    fit_params.add('T'+'_'+arg,value= tguess,vary=True,min=tguess - 10.0,max=tguess + 30.0)
    fit_params.add('L'+'_'+arg,value= 1e12,min=1e5,max=1e14)

  cov_ss_1d = minimize(simultaneous_stack_sed_oned, fit_params, 
    args=(cfits_flat,), kws={'data1d':flat_maps,'err1d':flat_noise,'wavelengths':cwavelengths,'LenLayers':LenLayers,'zed':zed})

  v = cov_ss_1d.params.valuesdict()

  beta = np.asarray(v['b'])
  for ised in range(nlists):
    arg = lists[ised]
    arg=arg.replace('.','p')
    arg=arg.replace('-','_')
    Temp = np.asarray(v['T'+'_'+arg])
    Lir = np.asarray(v['T'+'_'+arg])
    #Temp = np.asarray(v['T'+str(ised)])
    #Lir = np.asarray(v['L'+str(ised)])
    #pdb.set_trace()
    stacked_sed[:,ised]=single_simple_flux_from_greybody(np.sort(np.asarray(cwavelengths)), Trf = Temp, Lrf = Lir, b=beta, zin=zed)
  #pdb.set_trace()
  return [stacked_sed, v]
