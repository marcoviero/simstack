#!/usr/bin/env python

# Standard modules
import pdb
import os.path
import sys
import shutil
import time
import logging
import importlib
from astropy.wcs import WCS

# Modules within this package
import parameters 
from utils import circle_mask
from utils import dist_idl
from utils import gauss_kern
from utils import pad_and_smooth_psf
from utils import shift_twod
from utils import smooth_psf
from lmfit import Parameters, minimize, fit_report

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
    maps, psfs = get_maps(params)
    cats       = get_cats(params)
    bins       = get_bins(cats, params)

    # Do simultaneous stacking 

    fluxes	   = do_simstack(maps, psfs, cats, bins)

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

	return maps, psfs

if __name__=="__main__":
    main()
else:
    logging.info("Note: `mapit` module not being run as main executable.")
