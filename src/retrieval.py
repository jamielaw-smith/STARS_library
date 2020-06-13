"""
Aim:
One should be able to specify a (list of) parameters (mass, age, beta) to retrieve a dM/dt curve for.
I think we want to implement this as a quick 3D interpolation for each row in the list.
"""

import os
import numpy as np
from scipy.interpolate import interp1d
import shutil

def interpolate_beta(
    dmdt_input_dir = '../output/',
    dmdt_sub_dir = 'm0.925_t0.0/',
    output_dir = '../retrieval_scratch/',
    sim_beta_files = ['1.032.dat', '1.137.dat'], 
    beta_arr = np.array([1.1])):
    
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
    Sim_beta = [float(b[:-4]) for b in sim_beta_files]

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


def interpolate_mass(
    model_dir_formatter = "{}_{}",
    dmdt_input_dir = '../output/',
    #dmdt_sub_dir = 'm0.925_t0.0/',
    output_dir = '../retrieval_scratch/',
    #sim_beta_files = ['1.032.dat', '1.137.dat'],
    mass_steps = ['m0.3', 'm0.35'],
    age_strings = ['t0.0', 't1.0'],
    mass_arr = np.array([0.31]),
    #beta_arr = np.array([1.1])
    ):
    
    for age_string in age_strings:
        for m1, m2 in zip(mass_steps[:-1], mass_steps[1:]):

            print("Processing: [%s, %s] for time %s" % (m1, m2, age_string))

            m1_dir = model_dir_formatter.format(m1, age_string)
            m2_dir = model_dir_formatter.format(m2, age_string)
            dmdt_sub_dirs = [m1_dir, m2_dir]

            # --------- GET SIMULATION BETAS -----------------
            lo_sim_beta_files = os.listdir(dmdt_input_dir + dmdt_sub_dirs[0])
            lo_sim_beta_files.sort()
            lo_sim_betas = [float(b[:-4]) for b in lo_sim_beta_files]

            # find lo_interp_beta_files
            a = [f for f in os.listdir(dmdt_input_dir + dmdt_sub_dirs[0]) if not f.startswith('.')]
            a.sort()
            lo_interp_beta_files = a
            lo_interp_betas = [float(b[:-4]) for b in lo_interp_beta_files]

            # TODO note JLS I changed this from m1[1:4] to m1[1:] here because didn't work for m0.35
            Sim_mass = [float(m1[1:]), float(m2[1:])]

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
                time['lo'], dmdt['lo'] = np.genfromtxt(dmdt_input_dir + dmdt_sub_dirs[0] + "/" + low_interp_beta_file, skip_header=1, unpack=True)
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
                    a = [f for f in os.listdir(dmdt_input_dir + dmdt_sub_dirs[i]) if not f.startswith('.')]
                    a.sort()
                    # this will select 11 betas from lowest to highest. spaced by log from logsapce
                    # hi_interp_beta_files = a[::10]
                    # hi_interp_beta_files.append(a[-1])
                    hi_interp_beta_files = a
                    hi_interp_betas = [float(b[:-4]) for b in hi_interp_beta_files]


                    if len(hi_interp_beta_files) != len(lo_interp_beta_files):
                        print('ERROR not same number of betas. Not sure what to do?')

                    time['hi'], dmdt['hi'] = np.genfromtxt(dmdt_input_dir + dmdt_sub_dirs[i] + "/" + hi_interp_beta_files[z], skip_header=1, unpack=True)
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
                            #print(Sim_mass)
                        else:
                            slope = ((dmdt['hi'][j] - dmdtinterp) / (Sim_mass[i] - Sim_mass[i - 1]))
                            #print(Sim_mass)

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

                        file_to_save = output_dir + savesmalldir + beta_string + '.dat'
                        print("\tSaving %s" % file_to_save)
                        np.savetxt(file_to_save,
                                   np.transpose([np.concatenate([timeinterpolated[0], timeinterpolated[1]]),
                                np.concatenate([dmdtinterpolated[0], dmdtinterpolated[1]])]))


def interpolate_age(
    model_dir_formatter = "{}_{}",
    dmdt_input_dir = '../output/',
    output_dir = '../retrieval_scratch/',
    mass_string = 'm0.31',
    age_steps = ['t0.0', 't1.0'],
    age_arr = np.array([0.44]),
    ):

    # Iterate over all mass entries and do the age interpolation
    #for mass_string, age_steps in  model_directories_by_mass.items():
    for t1, t2 in zip(age_steps[:-1], age_steps[1:]):

        print("Processing: [%s, %s] for mass %s" % (t1, t2, mass_string))

        t1_dir = model_dir_formatter.format(mass_string, t1)
        t2_dir = model_dir_formatter.format(mass_string, t2)
        dmdt_sub_dirs = [t1_dir, t2_dir]


        # --------- GET SIMULATION BETAS -----------------

        # find lo_interp_beta_files
        a = [f for f in os.listdir(dmdt_input_dir + dmdt_sub_dirs[0]) if not f.startswith('.')]
        a.sort()
        # this will select 11 betas from lowest to highest. spaced by log from logsapce
        # lo_interp_beta_files = a[::10]
        # lo_interp_beta_files.append(a[-1])
        lo_interp_beta_files = a
        lo_interp_betas = [float(b[:-4]) for b in lo_interp_beta_files]

        Sim_age = [float(f.split('t')[1].split('/')[0]) for f in dmdt_sub_dirs]


        #for z, low_sim_beta_file in enumerate(lo_sim_beta_files):
        for z, low_interp_beta_file in enumerate(lo_interp_beta_files):

        #for x, dmdtsmalldir in enumerate(dmdtsmalldirs):
            #if x+1 == len(dmdtsmalldirs):
                # at max mass
                # probably already wrote this mass directory in interpolated_dmdts directory
            #    print('here1')
            #    break

            #dmdtdir = dmdtbigdir + dmdtsmalldir

            # ------ DIRECTORY PARAMETERS -------

            # It is assumed that there are different files for each beta
            # (such as 2.500.dat for beta = 2.5)
            # The first row is energy, the second is dmde.

            # dictionaries with gamma's as keys.
            Beta_slope = [] # {gammas[0]: [], gammas[1]: []}
            Beta_yinter = [] # {gammas[0]: [], gammas[1]: []}
            #Sim_beta = [] #{gammas[0]: [], gammas[1]: []}
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
            time['lo'], dmdt['lo'] = np.genfromtxt(dmdt_input_dir + dmdt_sub_dirs[0] + "/" + low_interp_beta_file, skip_header=1, unpack=True)
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
            for i in range(1, len(Sim_age)):
                #hi_sim_beta_files = os.listdir(dmdtbigdir + dmdtsmalldirs[i])
                #hi_sim_beta_files.sort()
                a = [f for f in os.listdir(dmdt_input_dir + dmdt_sub_dirs[i]) if not f.startswith('.')]
                a.sort()
                # this will select 11 betas from lowest to highest. spaced by log from logsapce
                # hi_interp_beta_files = a[::10]
                # hi_interp_beta_files.append(a[-1])
                hi_interp_beta_files = a
                hi_interp_betas = [float(b[:-4]) for b in hi_interp_beta_files]
                #hi_sim_betas = [float(b[:-4]) for b in hi_sim_beta_files]

                if len(hi_interp_beta_files) != len(lo_interp_beta_files):
                    print('ERROR not same number of betas. Not sure what to do?')


                #time['hi'], dmdt['hi'] = np.genfromtxt(dmdtdir + sim_beta_files[i], skip_header=1, unpack=True)
                #time['hi'], dmdt['hi'] = np.genfromtxt(dmdtbigdir + dmdtsmalldirs[x+1] + hi_sim_beta_files[i], skip_header=1, unpack=True)
                time['hi'], dmdt['hi'] = np.genfromtxt(dmdt_input_dir + dmdt_sub_dirs[i] + "/" + hi_interp_beta_files[z], skip_header=1, unpack=True)

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
                        # todo probably need to edit this for mass
                        #slope = ((dmdtinterp - dmdt['lo'][j]) / (Sim_beta[i] - Sim_beta[i - 1]))
                        slope = ((dmdtinterp - dmdt['lo'][j]) / (Sim_age[i] - Sim_age[i - 1]))
                    else:
                        slope = ((dmdt['hi'][j] - dmdtinterp) / (Sim_age[i] - Sim_age[i - 1]))

                    # this is really Mass_slope now?
                    Beta_slope[-1].append(slope)

                    yinter1 = (dmdt[nointerp][j] - Beta_slope[-1][j] * Sim_age[i - 1])
                    yinter2 = (dmdtinterp - Beta_slope[-1][j] * Sim_age[i])
                    Beta_yinter[-1].append((yinter1 + yinter2) / 2.0)
                    Mapped_time[-1].append(np.array(mapped_time[nointerp][j]))

                time['lo'] = np.copy(time['hi'])
                dmdt['lo'] = np.copy(dmdt['hi'])


            #interp_index_low = [0 for i in beta_arr]
            #interp_index_high = [0 for i in beta_arr]
            interp_index_low = [0 for i in age_arr]
            interp_index_high = [0 for i in age_arr]
            #for i, b in enumerate(beta_arr):#range(len(beta_arr)):
            for i, b in enumerate(age_arr):#range(len(beta_arr)):
                #for j in range(len(Sim_beta)):
                for j in range(len(Sim_age)):
                    if b == Sim_age[j]:
                        # no need to interp, already have dmdt & t for this beta
                        beta_interp = False
                        interp_index_high[i] = j
                        interp_index_low[i] = j
                        break

                    if b < Sim_age[j]:
                        interp_index_high[i] = j
                        interp_index_low[i] = j - 1
                        beta_interp = True
                        break

            #if not os.path.exists(savebigdir + dmdtsmalldir):
            #    os.makedirs(savebigdir + dmdtsmalldir)

            for i in range(len(age_arr)):
                if interp_index_high[i] == interp_index_low[i]:
                    # already have this mass
                    # probably already wrote this mass directory in interpolated_dmdts directory
                    pass
                    #time, dmdt = np.genfromtxt(dmdtdir + '{0:f}'.format(beta_arr[i])[:5]+ '.dat', skip_header=1, unpack=True)
                    #np.savetxt(savebigdir + dmdtsmalldir + '{0:f}'.format(beta_arr[i])[:5]+ '.dat',
                    #            np.transpose([time, dmdt]))
                else:

                    dmdtinterpolated = np.array([
                        Beta_yinter[interp_index_low[i]][0] +
                        Beta_slope[interp_index_low[i]][0] * age_arr[i],
                        Beta_yinter[interp_index_low[i]][1] +
                        Beta_slope[interp_index_low[i]][1] * age_arr[i]])

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
                            (age_arr[i] -
                             Sim_age[interp_index_low[i]]) /
                            (Sim_age[interp_index_high[i]] -
                             Sim_age[interp_index_low[i]]))

                    timeinterpolated = np.array(timeinterpolated)

                    #print(np.shape(dmdtinterpolated, dmdtinterpolated[0]))

                    savesmalldir = mass_string + '_t' + str(round(age_arr[i], 3))[:5] + '/'
                    if not os.path.exists(output_dir + savesmalldir):
                        os.makedirs(output_dir + savesmalldir)


                    # todo maybe need to assign beta string earlier above?
                    beta_float = round(lo_interp_betas[z] + ((age_arr[i] - Sim_age[0]) / (Sim_age[-1] - Sim_age[0]))\
                                  * (hi_interp_betas[z] - lo_interp_betas[z]),3)

                    beta_string = '{0:f}'.format(beta_float)[:5]
                    file_to_save = output_dir + savesmalldir + beta_string + '.dat'
                    print("\tSaving %s" % file_to_save)

                    np.savetxt(file_to_save,
                               np.transpose([np.concatenate([timeinterpolated[0], timeinterpolated[1]]),
                            np.concatenate([dmdtinterpolated[0], dmdtinterpolated[1]])]))

inputs = [
#mass,     age,    beta
#[0.31,    0.0,    0.99]
#[0.31,    0.44,    0.99]
[0.31,    0.44,    0.918]
]

#for row in inputs:

# just try for one first, then TODO make this loop through all rows in input
row = inputs[0]
mass = row[0]
age = row[1]
beta = row[2]

# todo have some check for whether any of the masses are below the lowest mass in interpolated grid. also if age is outside of [0,1]


# look if mass is already in output
outputdirs = os.listdir('../output/')
outputdirs = np.array(outputdirs)
outputdirs = outputdirs[np.where(outputdirs != '.DS_Store')]

# sort outputdirs based on mass. just outputdirs.sort() doesn't work
m_array = []
for i, outputdir in enumerate(outputdirs):
    m_array.append(float(outputdir.split('_')[0][1:]))

outputdirs = [x for _,x in sorted(zip(m_array,outputdirs))]
print(outputdirs)

# ^ looks like this is already sorted by time for each mass.

found_mass = False
found_age = False
found_beta = False

beta_lo = np.nan
beta_hi = np.nan
age_lo = np.nan
age_hi = np.nan
mass_lo = np.nan
mass_hi = np.nan

# find closest neighbors, check if we already have anything
for i, outputdir in enumerate(outputdirs):

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

            for k, betadir in enumerate(betadirs):
                b = float(betadir.split('.dat')[0])

                if (k==1) & (b > beta):
                    print('ERROR: Requested beta below lowest beta for this stellar mass and age.\
                         Lowest beta is:')
                    print(b)
                    exit()

                elif b == beta:
                    print('already have this beta')
                    found_beta = True

                    print('Your requested dM/dt is here:')
                    print('../output/' + outputdir + '/' + betadir)
                    break

                elif b > beta:
                    print('b > beta')
                    print(k-1, k)
                    print(betadirs[k-1], betadirs[k])
                    beta_lo = betadirs[k-1]
                    beta_hi = betadirs[k]
                    break

            break

        elif t > age:
            print('t > age')
            print(i-1, i)
            print(outputdirs[i-1], outputdirs[i])
            age_lo = float(outputdirs[i-1].split('_')[1][1:])
            age_hi = float(outputdirs[i].split('_')[1][1:])
            break

    elif m > mass:
        print('m > mass')
        print(i-1, i)
        print(outputdirs[i-1], outputdirs[i])
        mass_lo = float(outputdirs[i-1].split('_')[0][1:])
        mass_hi = float(outputdirs[i].split('_')[0][1:])
        break

print(mass_lo, mass_hi)
print(age_lo, age_hi)
print(beta_lo, beta_hi)


if (found_mass == True) & (found_age == True) & (found_beta == True):
    print('We are done.')
    exit()

if (found_mass == True) & (found_age == True) & (found_beta == False):
    print('Found mass, age, not beta, handled above.')
    exit()

if (found_mass == True) & (found_age == False) & (found_beta == False):
    print('Found mass, not age or beta.')
    ## first, at the given mass, interpolate in age between the 2 closest ages
    
    #interpolate_age()

    ## then interpolate in beta between 2 closest betas as above
    # (check if already have beta)


if (found_mass == False) & (found_age == False) & (found_beta == False):
    print('Found nothing.')
    ## first, interpolate in mass betwen the 2 closest masses. 
    # at what age(s)?
    # t0.0 and t1.0?
    # or the closest existing ages from first interpolation, like t0.2 and t0.3?
    print('INTERPOLATE MASS')
    interpolate_mass(
        model_dir_formatter = "{}_{}",
        dmdt_input_dir = '../output/',
        output_dir = '../retrieval_scratch/',
        #mass_steps = ['m0.3', 'm0.35'],
        mass_steps = ['m' + str(mass_lo), 'm' + str(mass_hi)],
        age_strings = ['t0.0', 't1.0'],
        mass_arr = np.array([mass]),
        )

    # check if already have age
    newdirs = os.listdir('../retrieval_scratch/')
    print(newdirs)
    if 'm' + str(mass) + '_' + 't' + str(age) not in newdirs:
        # don't already have age
        ## second, interpolate in age between the 2 closest ages as above
        print('INTERPOLATE AGE')
        interpolate_age(
            model_dir_formatter = "{}_{}",
            dmdt_input_dir = '../retrieval_scratch/',   # hmm TODO this should change depending if we have mass or not
            output_dir = '../retrieval_scratch/',
            mass_string = 'm' + str(mass),
            age_steps = ['t0.0', 't1.0'],
            age_arr = np.array([age]),
            )

    # check if already have beta
    betadirs = os.listdir('../retrieval_scratch/' + 'm' + str(mass) + '_' + 't' + str(age) + '/')
    if str(round(beta, 4)).ljust(5, '0') + '.dat' in betadirs:
        # copy existing beta to retrieval
        if not os.path.exists('../retrieval/' + 'm' + str(mass) + '_' + 't' + str(age)):
            os.makedirs('../retrieval/' + 'm' + str(mass) + '_' + 't' + str(age))

        shutil.copyfile('../retrieval_scratch/' + 'm' + str(mass) + '_' + 't' + str(age) + '/' + str(round(beta, 4)).ljust(5, '0') + '.dat',
                        '../retrieval/' + 'm' + str(mass) + '_' + 't' + str(age) + '/' + str(round(beta, 4)).ljust(5, '0') + '.dat')

    else:    
        # don't already have beta
        ## then interpoalte in beta between 2 closest betas as above
        print('INTERPOLATE BETA')
        interpolate_beta(
            dmdt_input_dir = '../retrieval_scratch/',
            #dmdt_sub_dir = outputdir + '/',
            dmdt_sub_dir = 'm' + str(mass) + '_' + 't' + str(age) + '/',
            output_dir = '../retrieval/',
            #sim_beta_files = [betadirs[i-1], betadirs[i]], 
            sim_beta_files = ['0.917.dat', '1.025.dat'],
            beta_arr = np.array([beta])
            )   

# delete scratch directory tree
shutil.rmtree('../retrieval_scratch/')







