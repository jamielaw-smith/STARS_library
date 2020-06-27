# STARS_library

STARS: stellar TDEs (tidal disruption events) with abundances and realistic structures

## description:

STARS_library is a grid of fallback rates to the black hole (dM/dt) from 3D hydrodynamical simulations of tidal disruption events (TDEs) using realistic stellar models. Possible use cases are:

(1) One can download the dmdts from Law-Smith+2020a here: LINK.

(2) One can create an interpolated library of dmdts from these simulations, to arbitrary spacing in stellar mass, stellar age, and impact parameter.

(3) One can retrieve dmdts for any particular values of (stellar mass, stellar age, impact parameter), either at the command line or via a list in a file.


## installation:

`git clone https://github.com/jamielaw-smith/STARS_library.git`

## setup:

`conda env create -f environment.yml`

`pip run requirements.txt`


## usage:

`python src/main.py`

This creates an interpolated library in `output/`. 
The input and output directories and the parameters of the interpolation are set in `STARS.config`:
input_dir sets the input directory, AKA the backbone of simulated models.
output_dir sets the output directory, AKA where the interpolated library is output.
`NUM_BETA_INTERP_POINTS` sets the total number of beta’s for each stellar mass and age.
`NUM_MASS_INTERP_POINTS` sets the number of stellar masses in between each input mass.
`NUM_AGE_INTERP_POINTS` set the number of stellar ages in between each input age.

`python src/main.py -r 1.0 0.0 2.0`

This retrieves a dmdt for a star with 
mass=1.0 [M_sun],
age=0.0 [fractional; 0 == ZAMS, 1.0 == TAMS], 
beta=2.0 [impact parameter].
This command line option is appropriate for retrieving a single interpolated model at a time.
Results will be placed in retrieval/ by default. This can be customized in `STARS.config`:
`retrieval_input_dir` sets the directory the retrieval code looks for its backbone, AKA the interpolated library.
`retrieval_scratch_dir` sets the temporary scratch directory the retrieval creates and then deletes in its execution.
`retrieval_output_dir` sets the directory the retrieval outputs to.

`python src/main.py -g`

This retrieves a grid of dmdts from `RETRIEVE.par`.
One can specify a list of stellar masses, age, and impact parameters for which to retrieve dmdts in `RETRIEVE.par`, either tab- or comma-separated. The name of this file can be customized in `STARS.config`:
`retrieval_grid_file` sets the parameter file read for -g.

## notes:

If you use the `-r` or `-g` options and the interpolated library in `output/` has not been initialized yet, the code will do this automatically.

Errors will be thrown if you request a dmdt from a stellar mass, stellar age, or impact parameter outside the range of the interpolated library. Stellar mass is in the range [0.3,3] M_sun, stellar age [0.0, 1.0] fractional MS age, and impact parameter has a particular range for each star. The allowed range for beta is a little trickier to guess a priori, but the code will tell you what the range for a given star is if it throws an error.

