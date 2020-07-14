# STARS_library

STARS: stellar TDEs (tidal disruption events) with abundances and realistic structures

## Description

STARS_library is a grid of fallback rates to the black hole (dM/dt) from 3D hydrodynamical simulations of tidal disruption events (TDEs) using realistic stellar models. Possible use cases are:

(1) One can download the original dM/dt's from Law-Smith+2020a here (note these also exist in the `input/` folder of this repository): [.zip file](LINK) (6 MB).

(2) One can create an interpolated library of dM/dt's from these simulations, to arbitrary spacing in stellar mass, stellar age, and impact parameter. 

Several pre-packaged interpolated grids with different spacings are also available:

- beta's=10, stellar masses=21, stellar ages=5: [.zip file](LINK) (35 MB).

- beta's=1000, stellar masses=6, stellar ages=2: [.zip file](LINK) (324 MB).

- beta's=100, stellar masses=56, stellar ages=12: [.zip file](LINK) (3 GB).

(3) One can retrieve dM/dt's for any particular values of (stellar mass, stellar age, impact parameter), either at the command line or via a list in a file.


## Installation

`git clone https://github.com/jamielaw-smith/STARS_library.git`

## Setup

`cd STARS_library`

To install the requirements in a new conda environment:

`conda env create -f environment.yml`

`conda activate stars_lib`

OR, to install requiements via pip:

`pip install -r requirements.txt`


## Usage

### Initialize interpolated library:

`python STARS_library.py`

This creates an interpolated library in `output/`. 
The input and output directories and the parameters of the interpolation are set in `STARS.config`:

`input_dir` sets the input directory, AKA the backbone of simulated models.

`output_dir` sets the output directory, AKA where the interpolated library is output.

`NUM_BETA_INTERP_POINTS` sets the total number of beta’s for each stellar mass and age.

`NUM_MASS_INTERP_POINTS` sets the number of stellar masses in between each input mass.

`NUM_AGE_INTERP_POINTS` set the number of stellar ages in between each input age.

### Retrieve a single dM/dt:

`python STARS_library.py -r 1.0 0.0 2.0`

This retrieves a dM/dt for a star with 

mass=1.0 [M_sun],

age=0.0 [fractional; 0 == ZAMS, 1.0 == TAMS], 

beta=2.0 [impact parameter].

This command line option is appropriate for retrieving a single interpolated model at a time. Results will be placed in `retrieval/` by default. This can be customized in `STARS.config`:

`retrieval_input_dir` sets the directory the retrieval code looks for its backbone, AKA the interpolated library.

`retrieval_scratch_dir` sets the temporary scratch directory the retrieval creates and then deletes in its execution.

`retrieval_output_dir` sets the directory the retrieval outputs to.

### Retrieve a list of dM/dt's:

`python STARS_library.py -g`

This retrieves a list of dM/dt's from `RETRIEVE.par`.

One can specify a list of stellar masses, age, and impact parameters for which to retrieve dM/dt's in `RETRIEVE.par`, either tab- or comma-separated. The name of this file can be customized in `STARS.config`:

`retrieval_grid_file` sets the parameter file read for `-g`.

## Notes

If you use the `-r` or `-g` options and the interpolated library in `output/` has not been initialized yet, the code will do this automatically.

Errors will be thrown if you request a dM/dt from a stellar mass, stellar age, or impact parameter outside the range of the interpolated library. Stellar mass is in the range [0.3,3] M_sun, stellar age [0.0, 1.0] fractional MS age, and impact parameter has a particular range for each star. The allowed range for beta is a little trickier to guess a priori, but the code will tell you what the range for a given star is if it throws an error.

If you are renaming the default output directory, make sure that output directory (`output/` by default) is completely empty on the first initialize.

The dM/dt's in the .dat files were constructed with 200 data points per log interval in time in order to keep the final size of the interpolated library small. Please contact us if you would like dM/dt files that are more finely spaced in time.

## Images

Volume rendering of a 1 Msun ZAMS star in a beta=1 encounter with a 10^6 Msun BH:

<img src="https://people.ucsc.edu/~lawsmith/theme/images/image1.png" width="800" alt="Simulation snapshot">

dM/dt's from STARS_library on a small grid:

<img src="https://people.ucsc.edu/~lawsmith/theme/images/mdots.png" width="800" alt="dM/dt's">

## Reference
Please cite our paper if you use this repository: LINK TO ADS.
