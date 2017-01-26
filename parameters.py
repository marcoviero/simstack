import pdb
import ConfigParser
import os
import logging
import pprint
import astropy.cosmology as ac
#from astropy.cosmology import Planck15 as cosmo

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
    raw_cosmo_params = dict(config.items('cosmology'))
    raw_io_params = dict(config.items('io'))
    raw_maps_to_stack_params = dict(config.items('maps_to_stack'))
    raw_map_path_params = dict(config.items('map_path'))
    raw_map_file_params = dict(config.items('map_file'))
    raw_noise_path_params = dict(config.items('map_path'))
    raw_noise_file_params = dict(config.items('noise_file'))
    raw_beams_params = dict(config.items('beams'))
    raw_color_correction_params = dict(config.items('color_correction'))
    raw_catalogs_params = dict(config.items('catalogs'))

    raw_io_params['param_file_path'] = os.path.abspath(param_file_path) # Store parameter file path

    params = get_general_params(raw_params) # Convert "raw" config dictionary to "organized" dictionary `params`
    params['io'] = get_io_parameters(raw_io_params)
    params['cosmo'] = get_cosmology_parameters(raw_cosmo_params)
    params['map_files'] = get_maps_parameters(raw_maps_to_stack_params,raw_map_path_params,raw_map_file_params)
    params['noise_files'] = get_maps_parameters(raw_maps_to_stack_params,raw_noise_path_params,raw_noise_file_params)
    params['psfs'] = get_beams_parameters(raw_maps_to_stack_params,raw_beams_params)
    #params['cats'] = get_catalogs_parameters(raw_catalogs_params)
    #params['bins'] = get_binning_parameters(raw_binning_params)

    logging.info("---------- PARAMETER VALUES ----------")
    logging.info("======================================")
    logging.info("\n" + pprint.pformat(params, indent=4) + "\n")

    pdb.set_trace()
    return params

def get_general_params(raw_params):
    params = {} # Initialize parameter dictionary

    # Type of galaxy split.  Default is UVJ star-forming / quiescent
    try:
        params['galaxy_splitting_scheme'] = raw_params['populations']
    except KeyError:
        params['galaxy_splitting_scheme'] = 'sf-qt'

    # If running bootstrap
    try:
        params['bootstrap'] = is_true(raw_params, 'bootstrap')
    except KeyError:
        params['bootstrap'] = False
    
    # Number of bootstraps
    if params['bootstrap'] == True:
        params['boots'] = float(raw_params['boots'])

    # If running bootstrap
    try:
        params['optimal_binning'] = is_true(raw_params, 'optimal_binning')
    except KeyError:
        params['optimal_binning'] = False
    
    # If running bootstrap
    try:
        params['bin_lookback_time'] = is_true(raw_params, 'bin_lookbackt')
    except KeyError:
        params['bin_lookback_time'] = False
    
    # If running bootstrap
    try:
        params['stack_all_z_at_once'] = is_true(raw_params, 'all_z_at_once')
    except KeyError:
        params['stack_all_z_at_once'] = False
    
    return params

def get_io_parameters(raw_params):
    io = {}

    io['output_folder']           = raw_params['output_folder']
    io['flux_densities_filename'] = raw_params['flux_densities_filename']
    io['boot_fluxes_filename']    = raw_params['boot_fluxes_filename']
    io['param_file_path']         = raw_params['param_file_path']

    return io

def get_cosmology_parameters(raw_params):
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

    return cosmo

def get_maps_parameters(raw_maps_to_stack_params,raw_map_path_params,raw_map_file_params):
    maps = {}

    for imap in raw_maps_to_stack_params:
        if bool(raw_maps_to_stack_params[imap].split()[1]) == True:
            maps[imap+''] = raw_map_path_params[imap] + raw_map_file_params[imap] 

    return maps

def get_beams_parameters(raw_maps_to_stack_params,raw_beams_params):
    psfs = {}

    for imap in raw_maps_to_stack_params:
        if bool(raw_maps_to_stack_params[imap].split()[1]) == True:
            if is_float(raw_beams_params[imap]) == True:
                psfs[imap+'_fwhm'] = float(raw_beams_params[imap])
            else:
                psfs[imap+'_beam_file'] = raw_beams_params[imap] 

            pdb.set_trace()
    return psfs

def get_color_correction_parameters(raw_color_correction_params):
    color_correction = {}

    return color_correction

def get_color_catalog_parameters(raw_catalog_params):
    catalog = {}

    return catalog

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

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

### FOR TESTING ###
if __name__=='__main__':
    import os, sys
    import pprint

    param_fp = sys.argv[1]
    print("")
    print("Testing %s on %s..." % (os.path.basename(__file__), param_fp))
    print("")
    pprint.pprint(get_params(param_fp))