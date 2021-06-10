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
# import cPickle as pickle
import pickle
from astropy.wcs import WCS

# Modules within this package
import parameters
from bincatalogs import Field_catalogs
from simstack import stack_libraries_in_layers
from simstack import stack_libraries_in_layers_w_background
from simstack import is_true
from bootstrap import Bootstrap

# Modules within Utils package
sys.path.append(os.path.join('..', 'Utils'))
from skymaps import Skymaps
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

    zkey = params['zkey']
    mkey = params['mkey']
    rkey = params['ra_key']
    dkey = params['dec_key']

    t0 = time.time()

    if params['bins']['bin_in_lookback_time'] == True:
        z_pref = 'lookt'
    else:
        z_pref = 'z'
    # Stack in Slices or ALL AT ONCE Choice made here
    if params['bins']['stack_all_z_at_once'] == True:
        n_slices = 1
    else:
        n_slices = len(params['bins']['z_nodes']) - 1

    # Save Parameter file in folder
    save_paramfile(params)

    for i in range(n_slices):
        if params['bins']['stack_all_z_at_once'] == True:
            j = None
            stacked_flux_density_key = 'all_' + z_pref
        else:
            j = i
            if params['bins']['bin_in_lookback_time'] == True:
                stacked_flux_density_key = '{:.2f}'.format(params['bins']['t_nodes'][j]) + '-' + '{:.2f}'.format(
                    params['bins']['t_nodes'][j + 1])
            else:
                stacked_flux_density_key = str(params['bins']['t_nodes'][j]) + '-' + str(
                    params['bins']['t_nodes'][j + 1])

        print(stacked_flux_density_key)
        # From parameter file read maps, psfs, cats, and divide them into bins
        sky_library = get_maps(params)
        cats = get_catalogs(params)
        if params['bootstrap'] == True:
            pcat = Bootstrap(cats.table)

        # Bootstrap Loop Starts here
        for iboot in np.arange(params['number_of_boots']) + params['boot0']:
            # stacked_flux_densities = {}
            if params['bootstrap'] == True:
                print('Running ' + str(int(iboot)) + ' of ' + str(int(params['boot0'])) + '-' + str(
                    int(params['boot0'] + params['number_of_boots'] - 1)) + ' bootstraps')

                pcat.perturb_catalog(perturb_z=params['perturb_z'])
                bootcat = Field_catalogs(pcat.pseudo_cat, zkey=zkey, mkey=mkey, rkey=rkey, dkey=dkey)
                binned_ra_dec = get_bin_radec(params, bootcat, single_slice=j)
                if params['save_bin_ids'] == False:
                    bin_ids = None
                else:
                    bin_ids = get_bin_ids(params, bootcat, single_slice=j)
                out_file_path = params['io']['output_folder'] + '/bootstrapped_fluxes/' + params['io']['shortname']
                out_file_suffix = '_' + stacked_flux_density_key + '_boot_' + str(int(iboot))
            else:
                binned_ra_dec = get_bin_radec(params, cats, single_slice=j)
                if params['save_bin_ids'] == False:
                    bin_ids = None
                else:
                    bin_ids = get_bin_ids(params, cats, single_slice=j)
                out_file_path = params['io']['output_folder'] + '/simstack_fluxes/' + params['io']['shortname']
                out_file_suffix = '_' + stacked_flux_density_key

            # Do simultaneous stacking
            if params['float_background'] == True:
                stacked_flux_densities = stack_libraries_in_layers_w_background(sky_library, binned_ra_dec)
            else:
                stacked_flux_densities = stack_libraries_in_layers(sky_library, binned_ra_dec)

            save_stacked_fluxes(stacked_flux_densities, params, out_file_path, out_file_suffix, IDs=bin_ids)
        # pdb.set_trace()

    # Summarize timing
    t1 = time.time()
    tpass = t1 - t0

    logging.info("Done!")
    logging.info("")
    logging.info("Total time                        : {:.4f} minutes\n".format(tpass / 60.))


def get_maps(params):
    '''
    Read maps and psfs and store into dictionaries
    '''
    sky_library = {}

    for t in params['library_keys']:
        sky = Skymaps(params['map_files'][t], params['noise_files'][t], params['psfs'][t + '_fwhm'],
                      color_correction=params['color_correction'][t], beam_area=params['psfs'][t + '_beam_area'])
        sky.add_wavelength(params['wavelength'][t])
        sky.add_fwhm(params['psfs'][t + '_fwhm'])
        sky_library[t] = sky
    return sky_library


def get_catalogs(params):
    # Formatting no longer needed as
    tbl = pd.read_table(params['catalogs']['catalog_path'] + params['catalogs']['catalog_file'], sep=',')

    tbl['ID'] = range(len(tbl))
    if 'sfg' in tbl.keys():
        pass
    elif 'CLASS' in tbl.keys():
        tbl['sfg'] = tbl['CLASS']

    zkey = params['zkey']
    mkey = params['mkey']
    rkey = params['ra_key']
    dkey = params['dec_key']
    catout = Field_catalogs(tbl, zkey=zkey, mkey=mkey, rkey=rkey, dkey=dkey)

    return catout


def get_bin_ids(params, cats, single_slice=None):
    if single_slice == None:
        z_nodes = params['bins']['z_nodes']
    else:
        z_nodes = params['bins']['z_nodes'][single_slice:single_slice + 2]
    m_nodes = params['bins']['m_nodes']

    if params['galaxy_splitting_scheme'] == 'sf-qt':
        cats.separate_sf_qt()
        cats.get_sf_qt_mass_redshift_bins(z_nodes, m_nodes)
        bin_ids = cats.id_z_ms
    elif params['galaxy_splitting_scheme'] == '5pops':
        Fcut = params['cuts']['fcut']
        MIPS24_cut = params['cuts']['mips24_cut']
        cats.separate_5pops(Fcut=Fcut, MIPS24_cut=MIPS24_cut)
        cats.get_5pops_mass_redshift_bins(z_nodes, m_nodes)
        bin_ids = cats.id_z_ms_5pop
    elif params['galaxy_splitting_scheme'] == '4pops':
        Fcut = params['cuts']['fcut']
        age_cut = params['cuts']['age_cut']
        cats.separate_4pops(Fcut=Fcut, age_cut=age_cut)
        cats.get_4pops_mass_redshift_bins(z_nodes, m_nodes)
        bin_ids = cats.id_z_ms_4pop
    elif params['galaxy_splitting_scheme'] == 'uvj':
        c_nodes = params['populations']['c_nodes']
        c_names = params['populations']['pop_names']
        cats.table['UVJ'] = np.sqrt((cats.table['rf_U_V'] - np.min(cats.table['rf_U_V'])) ** 2 + (
                cats.table['rf_V_J'] - np.min(cats.table['rf_V_J'])) ** 2)
        cats.separate_uvj_pops(c_nodes)
        cats.get_mass_redshift_uvj_bins(z_nodes, m_nodes, c_names)
        bin_ids = cats.id_z_ms_pop
    elif params['galaxy_splitting_scheme'] == 'general':
        cuts_dict = params['populations']
        cats.separate_pops_by_name(cuts_dict)
        cats.get_subpop_ids(z_nodes, m_nodes, cuts_dict)
        bin_ids = cats.subpop_ids

    return bin_ids


def get_bin_radec(params, cats, single_slice=None):
    if single_slice == None:
        z_nodes = params['bins']['z_nodes']
    else:
        z_nodes = params['bins']['z_nodes'][single_slice:single_slice + 2]
    m_nodes = params['bins']['m_nodes']

    if params['galaxy_splitting_scheme'] == 'sf-qt':
        cats.separate_sf_qt()
        cats.get_sf_qt_mass_redshift_bins(z_nodes, m_nodes)
        binned_ra_dec = cats.subset_positions(cats.id_z_ms)
    elif params['galaxy_splitting_scheme'] == '5pops':
        Fcut = params['cuts']['fcut']
        MIPS24_cut = params['cuts']['mips24_cut']
        cats.separate_5pops(Fcut=Fcut, MIPS24_cut=MIPS24_cut)
        cats.get_5pops_mass_redshift_bins(z_nodes, m_nodes)
        binned_ra_dec = cats.subset_positions(cats.id_z_ms_5pop)
    elif params['galaxy_splitting_scheme'] == '4pops':
        Fcut = params['cuts']['fcut']
        age_cut = params['cuts']['age_cut']
        cats.separate_4pops(Fcut=Fcut, age_cut=age_cut)
        cats.get_4pops_mass_redshift_bins(z_nodes, m_nodes)
        binned_ra_dec = cats.subset_positions(cats.id_z_ms_4pop)
    elif params['galaxy_splitting_scheme'] == 'uvj':
        c_nodes = params['populations']['c_nodes']
        c_names = params['populations']['pop_names']
        cats.table['UVJ'] = np.sqrt((cats.table['rf_U_V'] - np.min(cats.table['rf_U_V'])) ** 2 + (
                cats.table['rf_V_J'] - np.min(cats.table['rf_V_J'])) ** 2)
        cats.separate_uvj_pops(c_nodes)
        cats.get_mass_redshift_uvj_bins(z_nodes, m_nodes, c_names)
        binned_ra_dec = cats.subset_positions(cats.id_z_ms_pop)
    elif params['galaxy_splitting_scheme'] == 'general':
        cuts_dict = params['populations']
        cats.separate_pops_by_name(cuts_dict)
        cats.get_subpop_ids(z_nodes, m_nodes, cuts_dict)
        binned_ra_dec = cats.subset_positions(cats.subpop_ids)

    print(z_nodes)
    return binned_ra_dec


def save_stacked_fluxes(stacked_fluxes, params, out_file_path, out_file_suffix, IDs=None):
    fpath = "%s/%s_%s%s.p" % (
        out_file_path, params['io']['flux_densities_filename'], params['io']['shortname'], out_file_suffix)
    print('pickling to ' + fpath)
    if not os.path.exists(out_file_path): os.makedirs(out_file_path)

    if IDs == None:
        pickle.dump(stacked_fluxes, open(fpath, "wb"))  # , protocol=2 )
    else:
        pickle.dump([IDs, stacked_fluxes], open(fpath, "wb"))  # , protocol=2 )


def save_paramfile(params):
    fp_in = params['io']['param_file_path']
    if params['bootstrap'] == True:
        outdir = params['io']['output_folder'] + '/bootstrapped_fluxes/' + params['io']['shortname']
    else:
        outdir = params['io']['output_folder'] + '/simstack_fluxes/' + params['io']['shortname']
    print('writing parameter file to ' + outdir)
    if not os.path.exists(outdir): os.makedirs(outdir)
    fname = os.path.basename(fp_in)
    fp_out = os.path.join(outdir, fname)

    logging.info("Copying parameter file...")
    logging.info("  FROM : {}".format(fp_in))
    logging.info("    TO : {}".format(fp_out))
    logging.info("")

    shutil.copyfile(fp_in, fp_out)


if __name__ == "__main__":
    main()
else:
    logging.info("Note: `mapit` module not being run as main executable.")
