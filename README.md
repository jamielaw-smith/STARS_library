[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4062018.svg)](https://doi.org/10.5281/zenodo.4062018)

# STARS_library

STARS: stellar TDEs (tidal disruption events) with abundances and realistic structures.

## Description

STARS_library is an interpolated grid of fallback rates to the black hole (dM/dt) from 3D hydrodynamical simulations of tidal disruption events (TDEs) using realistic stellar models. Possible use cases are:

(1) One can download the original dM/dt's from Law-Smith+2020a here (note these also exist in the `input/` folder of this repository): [.zip file](https://www.dropbox.com/s/k14hfp88n5e4kk0/STARS_library_input.zip?dl=1) (2 MB).

(2) One can create an interpolated library of dM/dt's from these simulations, to arbitrary spacing in stellar mass, stellar age, and impact parameter (see below). 

Pre-packaged interpolated grids with different spacings are also available:

- betas=10, stellar masses=21, stellar ages=5: [.zip file](https://www.dropbox.com/s/xohdcp5tylazsrg/STARS_library_output_10_5_5.zip?dl=1) (26 MB).

- betas=1000, stellar masses=6, stellar ages=2: [.zip file](https://www.dropbox.com/s/wbanglobc1xu385/STARS_library_output_1000_2_2.zip?dl=1) (273 MB).

- betas=100, stellar masses=56, stellar ages=12: [.zip file](https://www.dropbox.com/s/uibocgirikcw11s/STARS_library_output_100_12_12.zip?dl=1) (1.6 GB).

(3) One can retrieve dM/dt's for any particular values of (stellar mass, stellar age, impact parameter), either at the command line or via a list in a file (see below).


## Installation

`git clone https://github.com/jamielaw-smith/STARS_library.git`

## Setup

The requirements are minimal (roughly python 3, scipy, numpy), so you may not need the step below.

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

- `input_dir` sets the input directory, AKA the backbone of simulated models.

- `output_dir` sets the output directory, AKA where the interpolated library is output.

- `NUM_BETA_INTERP_POINTS` sets the total number of beta’s for each stellar mass and age.

- `NUM_MASS_INTERP_POINTS` sets the number of stellar masses in between each input mass.

- `NUM_AGE_INTERP_POINTS` set the number of stellar ages in between each input age.

### Retrieve a single dM/dt:

`python STARS_library.py -r 1.0 0.0 2.0`

This retrieves a dM/dt for a star with 

mass=1.0 [M_sun],

age=0.0 [fractional; 0.0 = ZAMS, 1.0 = min(10 Gyr, TAMS)], 

beta=2.0 [impact parameter].

This command line option is appropriate for retrieving a single interpolated model at a time. Results will be placed in `retrieval/` by default. This can be customized in `STARS.config`:

- `retrieval_input_dir` sets the directory the retrieval code looks for its backbone, AKA the interpolated library.

- `retrieval_scratch_dir` sets the temporary scratch directory the retrieval creates and then deletes in its execution.

- `retrieval_output_dir` sets the directory the retrieval outputs to.

### Retrieve a list of dM/dt's:

`python STARS_library.py -g`

This retrieves a list of dM/dt's from `RETRIEVE.par`.

One can specify a list of stellar masses, age, and impact parameters for which to retrieve dM/dt's in `RETRIEVE.par`, either tab- or comma-separated. The name of this file can be customized in `STARS.config`:

- `retrieval_grid_file` sets the parameter file read for `-g`.

## Notes

The columns of the .dat files are t (day), dm/dt (M_sun/yr). Time is relative to first pericenter passage.

If you use the `-r` or `-g` options and the interpolated library in `output/` has not been initialized yet, the code will do this automatically.

Errors will be thrown if you request a dM/dt from a stellar mass, stellar age, or impact parameter outside the range of the interpolated library. Stellar mass is in the range [0.1, 10.0] M_sun, stellar age [0.0 = ZAMS, 1.0 = min(10 Gyr, TAMS)] fractional MS age, and impact parameter has a particular range for each star. The code will tell you what the beta range for a given star is if it throws an error. One can also determine the allowed betas by looking at Table 1 in the paper or the contents of the `input/` and `output/` directories.

As discussed in the paper, for M_star < 0.8 M_sun, the oldest stellar age we ran is at 10 Gyr rather than TAMS (which is older than the age of the universe). Thus the stellar age number (0.0 to 1.0) for these stars is the fraction of 10 Gyr rather than the fractional MS age. (So as indicated above, fractional MS age is in [0.0 = ZAMS, 1.0 = min(10 Gyr, TAMS)]).

Rare issue: if you are renaming the default output directory, make sure that output directory (`output/` by default) is completely empty on the first initialize.

The dM/dt's in the .dat files were constructed with 200 data points per log interval in time in order to keep the final size of the interpolated library small. Please contact us if you would like dM/dt files that are more finely spaced in time.

The library uses a $10^6 M_\sun$ black hole. One can scale the dM/dt's with BH mass for nonrelativistic disruptions ($r_p > 10 r_g$) as follows (Eq. 3,4 in paper): $dM/dt \propto M_{\rm BH}^{-1/2}$, $t \propto M_{\rm BH}^{1/2}$. 

One can also extend to other stellar masses by scaling dM/dt's from stars with similar density profiles. For example, one can extend below $0.1 M_\sun$ as these stars all have approximately $\gamma=5/3$ stellar structures. The scalings with stellar mass and radius are (Eq. 3,4 in paper):
$dM/dt \propto M_\star^{2} R_\star^{-3/2}$,
$t \propto M_\star^{-1} R_\star^{3/2}$.

## Reference
Please see the paper for much more detail on the simulations and results. Please cite both 

(1) the paper: https://ui.adsabs.harvard.edu/abs/2020ApJ...905..141L. Here is the BibTex entry:

```
@ARTICLE{2020ApJ...905..141L,
       author = {{Law-Smith}, Jamie A.~P. and {Coulter}, David A. and {Guillochon}, James and {Mockler}, Brenna and {Ramirez-Ruiz}, Enrico},
        title = "{Stellar Tidal Disruption Events with Abundances and Realistic Structures (STARS): Library of Fallback Rates}",
      journal = {\apj},
     keywords = {Black hole physics, Active galaxies, Galaxy nuclei, Gravitation, Hydrodynamics, Main sequence stars, Tidal disruption, 159, 17, 609, 661, 1963, 1000, 1696, Astrophysics - High Energy Astrophysical Phenomena},
         year = 2020,
        month = dec,
       volume = {905},
       number = {2},
          eid = {141},
        pages = {141},
          doi = {10.3847/1538-4357/abc489},
archivePrefix = {arXiv},
       eprint = {2007.10996},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2020ApJ...905..141L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

(2) this tool: https://ui.adsabs.harvard.edu/abs/2020zndo...4062018L. Here is the BibTex entry:

```
@MISC{2020zndo...4062018L,
       author = {{Law-Smith}, Jamie A.~P. and {Coulter}, David A. and {Mockler}, Brenna},
        title = "{jamielaw-smith/STARS\_library v1.0.5}",
         year = 2020,
        month = oct,
          eid = {10.5281/zenodo.4062018},
          doi = {10.5281/zenodo.4062018},
      version = {v1.0.5},
    publisher = {Zenodo},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2020zndo...4062018L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

if you use this repository. Thanks!

## Issues / requests
If you encounter any issues, or if you have suggestions for improvements that would be useful to you, please do not hesitate to let us know. You can reach Jamie Law-Smith at <lawsmith@ucsc.edu>. 
