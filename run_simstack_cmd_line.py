#!/usr/bin/env python

# Standard modules
import pdb
import os
import os.path
import sys
import shutil
import time
import logging
import importlib
import numpy as np
import pandas as pd
import cPickle as pickle
from astropy.wcs import WCS

# Modules within this package
import parameters
from skymaps import Skymaps
from bincatalogs import Field_catalogs
from utils import circle_mask
from utils import dist_idl
from utils import gauss_kern
from utils import pad_and_smooth_psf
from utils import shift_twod
from utils import smooth_psf
from lmfit import Parameters, minimize, fit_report
from simstack import stack_libraries_in_layers
from bootstrap import Bootstrap

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

    if params['bins']['bin_in_lookback_time'] == True:
        z_pref = 'lookt'
    else:
        z_pref = 'z'
    # Stack in Slices or ALL AT ONCE Choice made here
    if params['bins']['stack_all_z_at_once'] == True: n_slices = 1
    else: n_slices = len(params['bins']['z_nodes']) - 1

    for i in range(n_slices):
        if params['bins']['stack_all_z_at_once'] == True:
            j = None
            stacked_flux_density_key = 'all_'+z_pref
        else:
            j = i
            stacked_flux_density_key = str(params['bins']['t_nodes'][j])+'-'+str(params['bins']['t_nodes'][j+1])
        print stacked_flux_density_key
        # From parameter file read maps, psfs, cats, and divide them into bins
        sky_library   = get_maps(params)
        cats          = get_catalogs(params)
        if params['bootstrap'] == True:
            pcat = Bootstrap(cats.table)

        # Bootstrap Loop Starts here
        for iboot in np.arange(params['number_of_boots'])+params['boot0']:
            stacked_flux_densities = {}
            if params['bootstrap'] == True:
                print 'Running ' +str(int(iboot))+' of '+ str(int(params['boot0'])) +'-'+ str(int(params['boot0']+params['number_of_boots']-1)) + ' bootstraps'
                if params['index_boots'] == True:
                    boot_indices_path = params['io']['output_folder']+'/bootstrapped_fluxes/bootstrap_indices/'
                    boot_index_key = ''
                    if not os.path.exists(boot_indices_path): os.makedirs(boot_indices_path)
                else:
                    boot_indices_path = False
                    boot_index_key = None
                # pcat is an instance of Bootstrap
                # pcat.perturb_catalog with options encoded in parameter file
                pcat.perturb_catalog(perturb_z = params['perturb_z'], boot_indices_path = boot_indices_path, boot_index_key = boot_index_key)
                bootcat = Field_catalogs(pcat.pseudo_cat)
                #bootcat = Field_catalogs(Bootstrap(pcat.pseudo_cat).table)
                #bootcat = Field_catalogs(Bootstrap(cats.table).table)
                binned_ra_dec = get_bin_radec(params, bootcat, single_slice = j)
                #shortname = params['shortname']
                out_file_path   = params['io']['output_folder']+'/bootstrapped_fluxes/'+params['io']['shortname']
                out_file_suffix = '_'+stacked_flux_density_key+'_boot_'+str(int(iboot))
            else:
                binned_ra_dec = get_bin_radec(params, cats, single_slice = j)
                #shortname = params['shortname']
                out_file_path   = params['io']['output_folder'] + '/simstack_fluxes/' + params['io']['shortname']
                out_file_suffix = '_'+stacked_flux_density_key

            # Do simultaneous stacking
            #pdb.set_trace()
            stacked_flux_densities[stacked_flux_density_key] = stack_libraries_in_layers(sky_library,binned_ra_dec)

            save_stacked_fluxes(stacked_flux_densities,params,out_file_path,out_file_suffix)
        #pdb.set_trace()

    #Save Parameter file in folder
    save_paramfile(params)

    # Summarize timing
    t1 = time.time()
    tpass = t1-t0

    logging.info("Done!")
    logging.info("")
    logging.info("Total time                        : {:.4f} minutes\n".format(tpass/60.))

def get_maps(params):
    '''
    Read maps and psfs and store into dictionaries
    '''
    sky_library = {}

    for t in params['library_keys']:
        sky = Skymaps(params['map_files'][t],params['noise_files'][t],params['psfs'][t+'_fwhm'],color_correction=params['color_correction'][t], beam_area=params['psfs'][t+'_beam_area'])
        sky.add_wavelength(params['wavelength'][t])
        sky.add_fwhm(params['psfs'][t+'_fwhm'])
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
        #tbl['z_err']=tbl[['ZPDF_L68','ZPDF_H68']].mean(axis=1)
        tbl['z_err']=((tbl['ZPDF']-tbl['ZPDF_L68']) + (tbl['ZPDF_H68']-tbl['ZPDF']))/2.0

    if 'LMASS' in tbl.keys():
        pass
    elif 'MASS_MED' in tbl.keys():
        tbl['LMASS']=tbl['MASS_MED']
        #tbl['LMASS_ERR']=tbl[['MASS_MED_MIN68','MASS_MED_MAX68']].mean(axis=1)
        tbl['LMASS_ERR']=((tbl['MASS_MED']-tbl['MASS_MED_MIN68']) + (tbl['MASS_MED_MAX68']-tbl['MASS_MED']))/2.0

    if 'lmass' in tbl.keys():
        tbl['LMASS'] = tbl['lmass']

    if 'sfg' in tbl.keys():
        pass
    elif 'CLASS' in tbl.keys():
        tbl['sfg']=tbl['CLASS']

    catout = Field_catalogs(tbl)
    #try:
    #    catout.table['sfg']
    #    pass
    #except KeyError:
    #    if len(params['populations']) > 2:
    #        catout.separate_pops_by_name(params['populations'])
    #    else:
    #        catout.separate_sf_qt()

    return catout

def get_bin_radec(params, cats, single_slice = None):

    if single_slice == None:
        z_nodes = params['bins']['z_nodes']
    else:
        z_nodes = params['bins']['z_nodes'][single_slice:single_slice+2]
    m_nodes = params['bins']['m_nodes']

    #if params['populations'] > 2:
    #    cats.get_subpop_ids(z_nodes, m_nodes, params['populations'])
    #    binned_ra_dec = cats.subset_positions(cats.subpop_ids)
    #else:
    #    cats.get_sf_qt_mass_redshift_bins(z_nodes,m_nodes)
    #    binned_ra_dec = cats.subset_positions(cats.id_z_ms)

    if params['galaxy_splitting_scheme'] == 'sf-qt':
        cats.separate_sf_qt()
        cats.get_sf_qt_mass_redshift_bins(z_nodes,m_nodes)
        binned_ra_dec = cats.subset_positions(cats.id_z_ms)
    elif params['galaxy_splitting_scheme'] == '5pops':
        Fcut = params['cuts']['fcut']
        MIPS24_cut = params['cuts']['mips24_cut']
        cats.separate_5pops(Fcut=Fcut,MIPS24_cut=MIPS24_cut)
        cats.get_5pops_mass_redshift_bins(z_nodes,m_nodes)
        binned_ra_dec = cats.subset_positions(cats.id_z_ms_5pop)
    elif params['galaxy_splitting_scheme'] == '4pops':
        Fcut = params['cuts']['fcut']
        age_cut = params['cuts']['age_cut']
        cats.separate_4pops(Fcut=Fcut,age_cut=age_cut)
        cats.get_4pops_mass_redshift_bins(z_nodes,m_nodes)
        binned_ra_dec = cats.subset_positions(cats.id_z_ms_4pop)
    elif params['galaxy_splitting_scheme'] == 'uvj':
        c_nodes = params['populations']['c_nodes']
        c_names = params['populations']['pop_names']
        cats.table['UVJ']=np.sqrt((cats.table['rf_U_V'] - np.min(cats.table['rf_U_V']))**2 + (cats.table['rf_V_J']-np.min(cats.table['rf_V_J'])) ** 2)
        cats.separate_uvj_pops(c_nodes)
        cats.get_mass_redshift_uvj_bins(z_nodes,m_nodes,c_names)
        binned_ra_dec = cats.subset_positions(cats.id_z_ms_pop)
    elif params['galaxy_splitting_scheme'] == 'general':
        cuts_dict = params['populations']
        cats.separate_pops_by_name(cuts_dict)
        cats.get_subpop_ids(z_nodes, m_nodes, cuts_dict)
        binned_ra_dec = cats.subset_positions(cats.subpop_ids)

    print z_nodes
    return binned_ra_dec

def save_stacked_fluxes(stacked_fluxes, params, out_file_path, out_file_suffix):
    fpath = "%s/%s_%s%s.p" % (out_file_path, params['io']['flux_densities_filename'],params['io']['shortname'],out_file_suffix)
    print 'pickling to '+fpath
    if not os.path.exists(out_file_path): os.makedirs(out_file_path)

    nodes = params['bins']
    #pdb.set_trace()
    #np.savez(fpath, stacked_fluxes=stacked_fluxes, nodes=nodes)
    pickle.dump( [nodes, stacked_fluxes], open( fpath, "wb" ) )

def save_paramfile(params):
    fp_in    = params['io']['param_file_path']
    if params['bootstrap'] == True:
        outdir   = params['io']['output_folder']+'/bootstrapped_fluxes/'+params['io']['shortname']
    else:
        outdir   = params['io']['output_folder']+'/simstack_fluxes/'+params['io']['shortname']
    fname    = os.path.basename(fp_in)
    fp_out   = os.path.join(outdir, fname)

    logging.info("Copying parameter file...")
    logging.info("  FROM : {}".format(fp_in))
    logging.info("    TO : {}".format(fp_out))
    logging.info("")

    shutil.copyfile(fp_in, fp_out)

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
