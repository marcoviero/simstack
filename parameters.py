import ConfigParser
import os
import logging
import pprint
import astropy.cosmology as ac
from astropy.cosmology import Planck15 as cosmo

def get_params(param_file_path):
    """
    Get parameter values and return them in the form of a dictionary.

    Parameters
    ----------
    param_file_path : str
        path to parameter file

    Returns
    -------
    params : dict
    """

    config = ConfigParser.SafeConfigParser()
    config.read(param_file_path)

    # Get "raw" dictionaries from `config` object
    raw_params = dict(config.items('general'))
    raw_io_params = dict(config.items('io'))
    raw_cosmo_params = dict(config.items('cosmology'))
    raw_maps_params = dict(config.items('maps'))
    raw_beams_params = dict(config.items('beams'))
    raw_catalogs_params = dict(config.items('catalogs'))
    raw_binning_params = dict(config.items('binning'))
    raw_bootstrap_params = dict(config.items('bootstrap'))

    raw_io_params['param_file_path'] = os.path.abspath(param_file_path) # Store parameter file path

    params = get_general_params(raw_params) # Convert "raw" config dictionary to "organized" dictionary `params`
    params['io'] = get_io_parameters(raw_io_params)
    params['cosmo'] = get_cosmology_parameters(raw_cosmology_params)
    params['maps'] = get_maps_parameters(raw_maps_params)
    params['psfs'] = get_beams_parameters(raw_beams_params)
    params['cats'] = get_catalogs_parameters(raw_catalogs_params)
    params['bins'] = get_binning_parameters(raw_binning_params)
    params['boot'] = get_bootstrap_parameters(raw_bootstrap_params)

    logging.info("---------- PARAMETER VALUES ----------")
    logging.info("======================================")
    logging.info("\n" + pprint.pformat(params, indent=4) + "\n")

    return params

def get_io_parameters(raw_params):
    io = {}


    io['output_folder']         = raw_params['output_folder']
    io['fname_tcube']           = raw_params['fname_tcube']
    io['fname_powerspectrum']   = raw_params['fname_powerspectrum']
    try:
        io['save_tcube']        = is_true(raw_params, 'save_tcube')
    except KeyError:
        io['save_tcube']        = True

    io['param_file_path']       = raw_params['param_file_path']

    return io

def get_cosmology_parameters(raw_cosmology_params):
    '''
    Returns
    -------
    cosmo : astropy.cosmology object
        object containing cosmological parameters
    '''
    omega_m0    = float(raw_params['omega_m'])    # Present-day matter density
    omega_l0    = float(raw_params['omega_l'])    # Present-day dark energy density
    omega_k0    = float(raw_params['omega_k'])    # Present-day spatial curvature density
    hubble_h0   = float(raw_params['h'])          # Present-day reduced Hubble constant: h0 = H0 / (100 km/s/Mpc)

    H0          = hubble_h0*100.
    cosmo       = ac.LambdaCDM(H0=H0, Om0=omega_m0, Ode0=omega_l0)

    cosmo = {}

    return cosmo

def get_maps_parameters(raw_maps_params):
    maps = {}

    io['maps_path']        = raw_params['maps_path']
    return maps

def get_beams_parameters(raw_beams_params):
    psfs = {}

    return psfs

def get_catalogs_parameters(raw_catalogs_params):
    cats = {}

    return cats

def get_cosmology_parameters(raw_binning_params):
    bins = {}

    return bins

def get_cosmology_parameters(raw_bootstrap_params):
    boots = {}

    return boots
