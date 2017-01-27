#!/usr/bin/env python

# Standard modules
import pdb
import os.path
import sys
import shutil
import time
import logging
import importlib
import pandas as pd
from astropy.wcs import WCS

# Modules within this package
import parameters 
from skymaps import Skymaps
from skymaps import Field_catalogs
from utils import circle_mask
from utils import dist_idl
from utils import gauss_kern
from utils import pad_and_smooth_psf
from utils import shift_twod
from utils import smooth_psf
from lmfit import Parameters, minimize, fit_report
from simstack import stack_libraries_in_layers

def main():
	# Set up logging
    logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s  %(message)s',
            datefmt='%Y-%d-%m %I:%M:%S %p')

    # Get parameters from the provided parameter file
    param_file_path = sys.argv[1]
    params = parameters.get_params(param_file_path)

    t0 = time.time()

    # From parameter file read maps, psfs, cats, and divide them into bins
    sky_library   = get_maps(params)
    cats          = get_catalogs(params)
    binned_ra_dec = get_bins(params, cats)

    # Do simultaneous stacking 
    pdb.set_trace()
    stacked_flux_densities = stack_libraries_in_layers(sky_library,binned_ra_dec)

    # Summarize timing
    t1 = time.time()
    tpass = t1-t0

    logging.info("Done!")
    logging.info("")
    logging.info("Total time                        : {:.4f}\n".format(tpass))

def get_maps(params):
    '''
    Read maps and psfs and store into dictionaries
    '''
    sky_library = {}
    
    for t in params['library_keys']:
        sky = Skymaps(params['map_files'][t],params['noise_files'][t],params['psfs'][t+'_fwhm'],color_correction=params['color_correction'][t])
        sky.add_wavelength(params['wavelength'][t])
        sky.add_fwhm(params['psfs'][t+'_fwhm']) 
        #pdb.set_trace()
        sky_library[t] = sky 
    return sky_library 

def get_catalogs(params):

    # This formats the different catalogs to work with existing code.  
    # It's not the most elegant solution; better would be to make the code work with different column names.  
    tbl = pd.read_table(params['catalogs']['catalog_path']+params['catalogs']['catalog_file'],sep=',')

    if 'ID' in tbl.keys():
        pass
    elif 'ALPHA_J2000' in tbl.keys():
            tbl['ID']=range(len(tbl['ALPHA_J2000']))
    elif 'id' in tbl.keys():
            tbl['ID']=tbl['id']

    if 'ra' in tbl.keys():
        pass
    elif 'ALPHA_J2000' in tbl.keys():
        tbl['ra']=tbl['ALPHA_J2000']
        tbl['dec']=tbl['DELTA_J2000']

    if 'z_peak' in tbl.keys():
        pass
    elif 'ZPDF' in tbl.keys():
        tbl['z_peak']=tbl['ZPDF']

    if 'LMASS' in tbl.keys():
        pass
    elif 'MASS_MED' in tbl.keys():
        tbl['LMASS']=tbl['MASS_MED']

    if 'sfg' in tbl.keys():
        pass
    elif 'CLASS' in tbl.keys():
        tbl['sfg']=tbl['CLASS']

    return Field_catalogs(tbl)

def get_bins(params, cats):

    z_nodes = params['bins']['z_nodes']
    m_nodes = params['bins']['m_nodes']
    
    cats.get_sf_qt_mass_redshift_bins(z_nodes,m_nodes)
    binned_ra_dec = cats.subset_positions(cats.id_z_ms)

    return binned_ra_dec

def is_true(raw_params, key):
    """Is raw_params[key] true? Returns boolean value.
    """
    sraw    = raw_params[key]
    s       = sraw.lower() # Make case-insensitive

    # Lists of acceptable 'True' and 'False' strings
    true_strings    = ['true', 't', 'yes', 'y', '1']
    false_strings    = ['false', 'f', 'no', 'n', '0']
    if s in true_strings:
        return True
    elif s in false_strings:
        return False
    else:
        logging.warning("Input not recognized for parameter: %s" % (key))
        logging.warning("You provided: %s" % (sraw))
        raise

if __name__=="__main__":
    main()
else:
    logging.info("Note: `mapit` module not being run as main executable.")
