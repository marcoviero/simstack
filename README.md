# simstack
An open-source, Python version of SIMSTACK. 

For those familiar with the IDL version, viero_quick_stack is a similar wrapper to the stack_in_redshift_slices function.  
The python version of the code is intended to be object oriented, with map objects (containing maps, noisemaps, beams, color corrections, etc.) and catalog objects (with RA/DEC lists for different selection criteria) collected into libraries, and SIMSTACK performed from *the command line*, or from a script containing one line of code.  

## Getting Started

**As of *1/28/2017* SIMSTACK  got easier to use!** After setting local paths in your .bashrc (or .bash_profile), edit the parameter file (**use example.cfg as a guide**) to define your maps, catalogs, and pickled outputs (examples below), as well as to choose the binning and method of stacking,  and run from the command line with: 
../run_simstack_cmd_line.py example.cfg

There are just a few steps needed to get started. 
* Make sure all dependencies are installed
* Set your paths in the .bashrc (or .bash_profile)
* Edit the parameter file to 
	* contain paths and map/catalog names
	* desired binning
	* whether to stack in redshift bins, or all at once
	* map specific details like 
		* color corrections 
		* beam FWHM (solid angle if converting from MJy/sr to Jy/beam)

##### Dependencies
* Python 2
* Git repositories	
	* simstack
	* utils
	* [NDpredict](https://github.com/sawellons/NDpredict) *
	* [torrey_cmf](https://github.com/sawellons/NDpredict) *
* astropy
* numpy
* pandas
* [lmfit](https://lmfit.github.io/lmfit-py/index.html) (latest version required or it will fail!)

*_Optional repositories if you choose to bin masses by constant number densities._

##### Local Paths
In your .bashrc define the following, with edits to define your prefered directory locations:
 
	export MAPSPATH=$MAPSPATH/data/maps/
	export CATSPATH=$CATSPATH/data/catalogs/
	export PICKLESPATH=$PICKLESPATH/data/pickles/

## Setting up the Parameter file to run simstack from command line 
To-be-completed (but it's pretty self explainatory, follow example.cfg)
The only tricky part is defining [populations], which sets the the galaxy splitting scheme (and still needs work...)
The two recommended options (for now) are:
* sf-qt
	* Splits the sample into star-forming and quiescent according to the Williams UVJ definition.  Requires rest-frame U-J and V-J colors, labeled rf_U_V and rf_V_J, respectively.   
* general
	* Populations are defined by obeying __one__ condition (for now, ultimately would like no limit on conditions). To ignore condition set it to False.  
	Each line defines one population, and takes the form
		* population_key = population_index condition_key greater_that less_than equal_to
		For example: 
			+ sf = 1 CLASS False False 1
		the population 'sf' will be labeled 1 in the sfg column of the catalog, and is defined as objects with CLASS = 1 
        	+ sb = 2 lssfr -8.0 False
        	the population 'sb' will be labeled 2 in the sfg column of the catalog, and is defined as objects whose lssfr < -8.0 

	Conditions are set in reverse order, so if, say, a galaxy obeys the conditions of both sf and agn, but agn = 3 and sf = 1, then it will be defined as agn. 


## Code Outputs
Objects!  Tools to read and interpret those objects are under constant developement and improving all the time.  
[An iPython Notebook](https://github.com/marcoviero/simstack/tree/master/notebooks) is provided as an example of how to use the following functions:
* from simstack import PickledStacksReader
	* the workhorse.  Uses the info stored in the configuration file (which is rewritten to the output directory) to read in the stacks, determines if they are bootstraps, and if they are then calculates the uncertainties.  
* from simstack import measure_cib

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
