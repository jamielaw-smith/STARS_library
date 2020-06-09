"""
Aim:
One should be able to specify a (list of) parameters (mass, age, beta) to retrieve a dM/dt curve for.
I think we want to implement this as a quick 3D interpolation for each row in the list.
"""

import os
import numpy as np
from scipy.interpolate import interp1d






## try to define an interpolation function to call above
def interpolate_beta():
    print('here')
    print(outputdir)
    dmdt_sub_dir = 'm0.925_t0.0/'

    dmdt_input_dir = '../output/'
    output_dir = '../retrieval_scratch/'


    current_dmdt_dir = dmdt_input_dir + dmdt_sub_dir

    # dictionaries with gamma's as keys.
    Beta_slope = [] # {gammas[0]: [], gammas[1]: []}
    Beta_yinter = [] # {gammas[0]: [], gammas[1]: []}

    Mapped_time = [] # {gammas[0]: [], gammas[1]: []}
    # for converting back from mapped time to actual times and doing
    # interpolation in actual time
    Premaptime = [] #{gammas[0]: [], gammas[1]: []}
    Premapdmdt = [] #{gammas[0]: [], gammas[1]: []}


    # --------- GET SIMULATION BETAS -----------------

    #sim_beta_files = [betadirs[i-1], betadirs[i]]
    sim_beta_files = ['1.032.dat', '1.137.dat']

    Sim_beta = [float(b[:-4]) for b in sim_beta_files]

    beta_arr = np.array([1.1])

    # ----- CREATE INTERPOLATION FUNCTIONS; FIND SLOPES & YINTERs -----
    time = {}
    dmdt = {}
    ipeak = {}
    mapped_time = {}

    # get dmdt and t for the lowest beta value
    # energy & dmde (cgs)
    time['lo'], dmdt['lo'] = np.genfromtxt(current_dmdt_dir + sim_beta_files[0], skip_header=1, unpack=True)
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

    for i in range(1, len(Sim_beta)):
        # indexing this way bc calculating slope and yintercepts
        # BETWEEN each simulation beta

        time['hi'], dmdt['hi'] = np.genfromtxt(current_dmdt_dir + sim_beta_files[i], skip_header=1, unpack=True)

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
                slope = ((dmdtinterp - dmdt['lo'][j]) /
                         (Sim_beta[i] - Sim_beta[
                             i - 1]))
            else:
                slope = ((dmdt['hi'][j] - dmdtinterp) /
                         (Sim_beta[i] - Sim_beta[
                             i - 1]))
            Beta_slope[-1].append(slope)

            yinter1 = (dmdt[nointerp][j] - Beta_slope[-1][j] *
                       Sim_beta[i - 1])
            yinter2 = (dmdtinterp - Beta_slope[-1][j] *
                       Sim_beta[i])
            Beta_yinter[-1].append((yinter1 + yinter2) / 2.0)
            Mapped_time[-1].append(
                np.array(mapped_time[nointerp][j]))

        time['lo'] = np.copy(time['hi'])
        dmdt['lo'] = np.copy(dmdt['hi'])


    interp_index_low = [0 for i in beta_arr]
    interp_index_high = [0 for i in beta_arr]
    for i, b in enumerate(beta_arr):#range(len(beta_arr)):
        for j in range(len(Sim_beta)):
            if b == Sim_beta[j]:
                # no need to interp, already have dmdt & t for this beta
                beta_interp = False
                interp_index_high[i] = j
                interp_index_low[i] = j
                break

            if b < Sim_beta[j]:
                interp_index_high[i] = j
                interp_index_low[i] = j - 1
                beta_interp = True
                break

    if not os.path.exists(output_dir + dmdt_sub_dir):
        os.makedirs(output_dir + dmdt_sub_dir)

    for i in range(len(beta_arr)):
        if interp_index_high[i] == interp_index_low[i]:
            time, dmdt = np.genfromtxt(current_dmdt_dir + '{0:f}'.format(beta_arr[i])[:5] + '.dat', skip_header=1, unpack=True)
            np.savetxt(output_dir + dmdt_sub_dir + '{0:f}'.format(beta_arr[i])[:5] + '.dat',
                       np.transpose([time, dmdt]))
        else:

            dmdtinterpolated = np.array([
                                Beta_yinter[interp_index_low[i]][0] +
                                Beta_slope[interp_index_low[i]][0] * beta_arr[i],
                                Beta_yinter[interp_index_low[i]][1] +
                                Beta_slope[interp_index_low[i]][1] * beta_arr[i]])

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
                    (beta_arr[i] -
                     Sim_beta[interp_index_low[i]]) /
                    (Sim_beta[interp_index_high[i]] -
                     Sim_beta[interp_index_low[i]]))

            timeinterpolated = np.array(timeinterpolated)

            #print(np.shape(dmdtinterpolated, dmdtinterpolated[0]))
            np.savetxt(output_dir + dmdt_sub_dir + '{0:f}'.format(beta_arr[i])[:5] + '.dat', np.transpose([np.concatenate([timeinterpolated[0], timeinterpolated[1]]),
                                                                                                           np.concatenate([dmdtinterpolated[0], dmdtinterpolated[1]])]))








inputs = [
#mass,     age,    beta
[0.925,    0.0,    1.1]
]

#for row in inputs:

# just try for one first, then TODO make this loop through all rows in input
row = inputs[0]
mass = row[0]
age = row[1]
beta = row[2]

# look if mass is already in output
outputdirs = os.listdir('../output/')
outputdirs.sort()

found_mass = False
found_age = False
found_beta = False

for outputdir in outputdirs:
    if outputdir == '.DS_Store':
        continue

    m = float(outputdir.split('_')[0][1:])
    t = float(outputdir.split('_')[1][1:])

    if m == mass:
        print('already have this mass')
        found_mass = True

        if t == age:
            print('already have this age')
            found_age = True

            betadirs = os.listdir('../output/'+outputdir)
            betadirs.sort()
            for betadir in betadirs:
                b = float(betadir.split('.dat')[0])
                if b == beta:
                    print('already have this beta')
                    found_beta = True

                    print('Your requested dM/dt is here:')
                    print('../output/' + outputdir + '/' + betadir)

                    break

            if found_beta == False:
                for i, betadir in enumerate(betadirs):
                    b = float(betadir.split('.dat')[0])
                    if b < beta:
                        print('b < beta')
                    elif b > beta:
                        print('b > beta')
                        print(i-1, i)
                        print(betadirs[i-1], betadirs[i])
                        interpolate_beta()
                        break




if (found_mass == True) & (found_age == True) & (found_beta == True):
    print('We are done.')
    exit()

if (found_mass == True) & (found_age == True) & (found_beta == False):
    print('Found mass, age, not beta.')
    ## interpolate in beta between 2 closest betas
    
    #find closest 2 betas

    #call dmdt_interpolation_beta.py
    
    #create new outputs in ../retrieval_scratch/


if (found_mass == True) & (found_age == False) & (found_beta == False):
    print('Found mass, not age or beta.')
    ## first, at the given mass, interpolate in age between the 2 closest ages
    ## then interpolate in beta between 2 closest betas as above
    # (check if already have beta)


if (found_mass == False) & (found_age == False) & (found_beta == False):
    print('Found nothing.')
    ## first, interpolate in mass betwen the 2 closest masses

    ## second, interpolate in age between the 2 closest ages as above
    # (check if already have age)

    ## then interpoalte in beta between 2 closest betas as above
    # (check if already have beta)




## return a file(s) in new directory with requested parameters.








