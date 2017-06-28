# simstack
An open-source, Python version of SIMSTACK. 

**As of *1/28/2017* SIMSTACK  got easier to use!** Edit the parameter file (**use example.cfg as a guide**) to set the local paths to your maps, catalogs, and pickles (examples below), as well as to choose the binning and method of stacking  and run from the command line with: 
../run_simstack_cmd_line.py example.cfg

For those familiar with the IDL version, viero_quick_stack is a similar wrapper to the stack_in_redshift_slices function.  
The python version of the code is intended to be object oriented, with map objects (containing maps, noisemaps, beams, color corrections, etc.) and catalog objects (with RA/DEC lists for different selection criteria) collected into libraries, and SIMSTACK performed with one line of code.  

## Getting Started
##### Dependencies
* Git repositories	
	* simstack
	* utils
	* [NDpredict](https://github.com/sawellons/NDpredict) *
	* [torrey_cmf](https://github.com/sawellons/NDpredict) *
* astropy
* numpy
* pandas

_Optional repositories if you choose to bin masses by constant number densities._

##### Local Paths
In your .bashrc define the following, with edits to define your prefered directory locations:
 
	export MAPSPATH=$MAPSPATH/data/maps/
	export CATSPATH=$CATSPATH/data/catalogs/
	export PICKLESPATH=$PICKLESPATH/data/pickles/

##### Running simstack from Script 


## Future Work

Future projects include adding priors from SIMSTACK, and moving from a least-squares solution to marginalized probabilities using a full MCMC.

SEDSTACK is an attempt to fit entire SEDs to the full set of maps --- rather than fitting for flux densities one wavelength at a time --- in order to dig still deeper into the confusion-dominated data.  Outstanding issues are dealing with systematic offsets, covariance between maps, careful consideration of beam and pixel size differences and how the affect the chi-squared, and using actual beams instead of Gaussian approximations.  

Hopefully, with community input (which I encourage!) the code will grow in utility.

## People

* Marco Viero
* Lorenzo Moncelsi

## License, Credits etc

This is work in progress! If you use any of the code or ideas here in your research, please cite us as (Viero et al. in preparation).

All content Copyright 2015 the authors. The code is available for use under the MIT license.
