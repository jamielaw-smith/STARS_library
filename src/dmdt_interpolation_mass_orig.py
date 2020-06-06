"""Initialize module.
copied from _1param.py
interpolating in mass, mapping betas in between masses, keeping age constant
replaced _beta with _mass often
todo think about what to do for mdots past critical beta
todo right now it looks like works best when only run between two masses at a time
"""

import os
import astropy.constants as c
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import time as tm
start_time = tm.time()


C_CGS = c.c.cgs.value
FOUR_PI = 4*np.pi
M_SUN_CGS = c.M_sun.cgs.value

G = c.G.cgs.value  # 6.67259e-8 cm3 g-1 s-2
Mhbase = 1.0e6 * M_SUN_CGS  # this is the generic size of bh used

# USER INPUTS
NUM_MASS_INTERP_POINTS = 10

age_string = 't0.0' #'t1.0'

#only between two masses at a time for now, for a given age
dmdt_sub_dirs = [
# 'm0.3_t0.0/',
#'m0.3_t1.0/',
# 'm0.5_t0.0/',
#'m0.5_t1.0/',
#'m0.7_t0.0/',
#'m0.7_t1.0/',
'm1.0_t0.0/',
#'m1.0_t1.0/',
'm3.0_t0.0/',
#'m3.0_t1.0/',
]

dmdt_input_dir = '../input/'
output_dir = '../output/'


# --------- GET SIMULATION BETAS -----------------
lo_sim_beta_files = os.listdir(dmdt_input_dir + dmdt_sub_dirs[0])
lo_sim_beta_files.sort()
lo_sim_betas = [float(b[:-4]) for b in lo_sim_beta_files]

# find lo_interp_beta_files
a = [f for f in os.listdir(output_dir + dmdt_sub_dirs[0]) if not f.startswith('.')]
a.sort()

# this will select 11 betas from lowest to highest. spaced by log from logsapce
lo_interp_beta_files = a[::10]
lo_interp_beta_files.append(a[-1])
lo_interp_betas = [float(b[:-4]) for b in lo_interp_beta_files]


Sim_mass = [float(f[1:4]) for f in dmdt_sub_dirs]
# mass_arr = np.logspace(np.log10(Sim_mass[0]), np.log10(Sim_mass[-1]), num=NUM_MASS_INTERP_POINTS)
mass_arr = np.linspace(Sim_mass[0], Sim_mass[-1], num=NUM_MASS_INTERP_POINTS)

for z, low_interp_beta_file in enumerate(lo_interp_beta_files):
    # ------ DIRECTORY PARAMETERS -------

    # It is assumed that there are different files for each beta
    # (such as 2.500.dat for beta = 2.5)
    # The first row is energy, the second is dmde.

    # dictionaries with gamma's as keys.
    Beta_slope = [] # {gammas[0]: [], gammas[1]: []}
    Beta_yinter = [] # {gammas[0]: [], gammas[1]: []}

    Mapped_time = [] # {gammas[0]: [], gammas[1]: []}
    # for converting back from mapped time to actual times and doing
    # interpolation in actual time
    Premaptime = [] #{gammas[0]: [], gammas[1]: []}
    Premapdmdt = [] #{gammas[0]: [], gammas[1]: []}


    # ----- CREATE INTERPOLATION FUNCTIONS; FIND SLOPES & YINTERs -----
    time = {}
    dmdt = {}
    ipeak = {}
    mapped_time = {}
    # get dmdt and t for the lowest beta value
    # energy & dmde (cgs)
    time['lo'], dmdt['lo'] = np.genfromtxt(output_dir + dmdt_sub_dirs[0] + low_interp_beta_file, skip_header=1, unpack=True)
    ipeak['lo'] = np.argmax(dmdt['lo'])

    # split time['lo'] & dmdt['lo'] into pre-peak and post-peak array
    time['lo'] = np.array([
        time['lo'][:ipeak['lo']],
        time['lo'][ipeak['lo']:]])  # peak in array 2
    dmdt['lo'] = np.array([
        dmdt['lo'][:ipeak['lo']],
        dmdt['lo'][ipeak['lo']:]])  # peak in array 2

    # will contain time/dmdt arrays
    # (split into pre & post peak times/dmdts)
    # for each beta value
    Premaptime.append(np.copy(time['lo']))
    Premapdmdt.append(np.copy(dmdt['lo']))

    #for i in range(1, len(Sim_beta)):
    for i in range(1, len(Sim_mass)):
        a = [f for f in os.listdir(output_dir + dmdt_sub_dirs[i]) if not f.startswith('.')]
        a.sort()
        # this will select 11 betas from lowest to highest. spaced by log from logsapce
        hi_interp_beta_files = a[::10]
        hi_interp_beta_files.append(a[-1])
        hi_interp_betas = [float(b[:-4]) for b in hi_interp_beta_files]


        if len(hi_interp_beta_files) != len(lo_interp_beta_files):
            print('ERROR not same number of betas. Not sure what to do?')

        time['hi'], dmdt['hi'] = np.genfromtxt(output_dir + dmdt_sub_dirs[i] + hi_interp_beta_files[z], skip_header=1, unpack=True)
        ipeak['hi'] = np.argmax(dmdt['hi'])

        # split time_hi and dmdt_hi into pre-peak and post-peak array
        # peak in 2nd array
        time['hi'] = np.array([time['hi'][:ipeak['hi']],
                               time['hi'][ipeak['hi']:]])
        dmdt['hi'] = np.array([dmdt['hi'][:ipeak['hi']],
                               dmdt['hi'][ipeak['hi']:]])
        # will contain time/dmdt arrays
        # (split into pre & post peak times/dmdts)
        # for each beta value
        Premapdmdt.append(np.copy(dmdt['hi']))
        Premaptime.append(np.copy(time['hi']))

        mapped_time['hi'] = []
        mapped_time['lo'] = []

        Beta_slope.append([])

        Beta_yinter.append([])
        Mapped_time.append([])
        for j in [0, 1]:  # once before peak, once after peak
            # choose more densely sampled curve to map times to 0-1
            # less densely sampled curve will be interpolated to match
            if len(time['lo'][j]) < len(time['hi'][j]):
                # hi array more densely sampled
                interp = 'lo'
                nointerp = 'hi'
            else:
                # will also catch case where they have the same lengths
                interp = 'hi'
                nointerp = 'lo'
            # map times from more densely sampled curves
            # (both pre & post peak, might be from diff. dmdts)
            # to 0 - 1
            mapped_time[nointerp].append(
                1. / (time[nointerp][j][-1] - time[nointerp][j][0]) *
                (time[nointerp][j] - time[nointerp][j][0]))
            mapped_time[interp].append(
                1. / (time[interp][j][-1] - time[interp][j][0]) *
                (time[interp][j] - time[interp][j][0]))

            # ensure bounds are same for interp and nointerp
            # before interpolation
            # (should be 0 and 1 from above, but could be slightly off
            # due to rounding errors in python)
            mapped_time[interp][j][0] = 0
            mapped_time[interp][j][-1] = 1
            mapped_time[nointerp][j][0] = 0
            mapped_time[nointerp][j][-1] = 1

            func = interp1d(mapped_time[interp][j], dmdt[interp][j])
            dmdtinterp = func(mapped_time[nointerp][j])

            if interp == 'hi':
                slope = ((dmdtinterp - dmdt['lo'][j]) / (Sim_mass[i] - Sim_mass[i - 1]))
            else:
                slope = ((dmdt['hi'][j] - dmdtinterp) / (Sim_mass[i] - Sim_mass[i - 1]))

            # this is really Mass_slope now?
            Beta_slope[-1].append(slope)

            yinter1 = (dmdt[nointerp][j] - Beta_slope[-1][j] * Sim_mass[i - 1])
            yinter2 = (dmdtinterp - Beta_slope[-1][j] * Sim_mass[i])
            Beta_yinter[-1].append((yinter1 + yinter2) / 2.0)
            Mapped_time[-1].append(np.array(mapped_time[nointerp][j]))

        time['lo'] = np.copy(time['hi'])
        dmdt['lo'] = np.copy(dmdt['hi'])


    interp_index_low = [0 for i in mass_arr]
    interp_index_high = [0 for i in mass_arr]

    for i, b in enumerate(mass_arr):#range(len(beta_arr)):
        for j in range(len(Sim_mass)):
            if b == Sim_mass[j]:
                # no need to interp, already have dmdt & t for this beta
                beta_interp = False
                interp_index_high[i] = j
                interp_index_low[i] = j
                break

            if b < Sim_mass[j]:
                interp_index_high[i] = j
                interp_index_low[i] = j - 1
                beta_interp = True
                break

    for i in range(len(mass_arr)):
        if interp_index_high[i] == interp_index_low[i]:
            # already have this mass
            # probably already wrote this mass directory in interpolated_dmdts directory
            pass
        else:
            dmdtinterpolated = np.array([
                                Beta_yinter[interp_index_low[i]][0] +
                                Beta_slope[interp_index_low[i]][0] * mass_arr[i],
                                Beta_yinter[interp_index_low[i]][1] +
                                Beta_slope[interp_index_low[i]][1] * mass_arr[i]])

                            # map mapped_times back to actual times, requires interpolation
                            # in time
                            # first for pre peak times

            timeinterpolated = []
            for j in [0, 1]:
                # interp_index_low indexes beta
                # mapped time between beta low and beta high
                time_betalo = (
                    Mapped_time[interp_index_low[i]][j] *
                    (Premaptime[interp_index_low[i]][j][-1] -
                     Premaptime[interp_index_low[i]][j][0]) +
                    Premaptime[interp_index_low[i]][j][0])
                time_betahi = (
                    Mapped_time[interp_index_low[i]][j] *
                    (Premaptime[interp_index_high[i]][j][-1] -
                     Premaptime[interp_index_high[i]][j][0]) +
                    Premaptime[interp_index_high[i]][j][0])

                timeinterpolated.append(
                    time_betalo + (time_betahi - time_betalo) *
                    (mass_arr[i] -
                     Sim_mass[interp_index_low[i]]) /
                    (Sim_mass[interp_index_high[i]] -
                     Sim_mass[interp_index_low[i]]))

            timeinterpolated = np.array(timeinterpolated)

            savesmalldir = 'm' + str(round(mass_arr[i], 3))[:5] + '_' + age_string + '/'
            if not os.path.exists(output_dir + savesmalldir):
                os.makedirs(output_dir + savesmalldir)


            # todo maybe need to assign beta string earlier above?
            beta_float = round(lo_interp_betas[z] + ((mass_arr[i] - Sim_mass[0])/(Sim_mass[-1] - Sim_mass[0]))\
                          * (hi_interp_betas[z] - lo_interp_betas[z]),3)

            beta_string = '{0:f}'.format(beta_float)[:5]
            np.savetxt(output_dir + savesmalldir + beta_string + '.dat',
                       np.transpose([np.concatenate([timeinterpolated[0], timeinterpolated[1]]),
                    np.concatenate([dmdtinterpolated[0], dmdtinterpolated[1]])]))

stop_time = tm.time()
duration = stop_time - start_time
print("\n\nExecution: %s [sec]" % duration)