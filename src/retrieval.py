import os
import numpy as np
import shutil
from src.interpolation import beta_interpolate, mass_interpolate, age_interpolate

def retrieval(mass, age, beta, retrieval_input_dir, retrieval_scratch_dir, retrieval_output_dir):
    # check if mass is outside [0.3, 3]. also if age is outside of [0,1]
    if not (0.3 <= mass <= 3.0):
        raise Exception('ERROR: mass must be in [0.3, 3]')
    if not (0.0 <= age <= 1.0):
        raise Exception('ERROR: age must be in [0, 1]')

    # check if already have this in ../retrieval/, 
    if os.path.exists(retrieval_output_dir):
        dirs = [d for d in os.listdir(retrieval_output_dir) if not d.startswith('.')]
        m_array_temp = [float(d.split('_')[0][1:]) for d in dirs]
        t_array_temp = [float(d.split('_')[1][1:]) for d in dirs]

        for i, (m, t) in enumerate(zip(m_array_temp, t_array_temp)):
            
            if (mass == m) & (age == t):
                beta_files = [f for f in os.listdir(retrieval_output_dir + dirs[i]) if not f.startswith('.')]
                b_array = [float(f.split('.dat')[0]) for f in beta_files]
                beta_files = np.array(beta_files)
                b_array = np.array(b_array)
                sel_beta_temp = np.where(b_array == beta)[0]
                if len(sel_beta_temp) != 0:
                    print('Already retrieved %s' % retrieval_output_dir + dirs[i] + '/' + beta_files[sel_beta_temp][0])
                    return

    outputdirs = [d for d in os.listdir(retrieval_input_dir) if not d.startswith('.')]
    m_array = [float(d.split('_')[0][1:]) for d in outputdirs]
    t_array = [float(d.split('_')[1][1:]) for d in outputdirs]

    outputdirs = np.array(outputdirs)
    m_array = np.array(m_array)
    t_array = np.array(t_array)

    sel_mass_and_age = np.where((m_array == mass) & (t_array == age))[0]
    sel_mass = np.where(m_array == mass)[0]

    # check if already have this mass and age
    if len(sel_mass_and_age) != 0:
        print('#1')
        print('have mass and age:')
        print(outputdirs[sel_mass_and_age][0])

        # get betas for this dir
        beta_files = [f for f in os.listdir(retrieval_input_dir + outputdirs[sel_mass_and_age][0]) if not f.startswith('.')]
        b_array = [float(f.split('.dat')[0]) for f in beta_files]

        beta_files = np.array(beta_files)
        b_array = np.array(b_array)

        # check if beta is outside range for given directory
        if not (min(b_array) <= beta <= max(b_array)):
            print('ERROR: beta outside range for this mass and age, which is', [min(b_array),max(b_array)])
            exit()

        # check if already have beta
        sel_beta = np.where(b_array == beta)[0]
        if len(sel_beta) != 0:
            # copy existing beta to retrieval
            if not os.path.exists(retrieval_output_dir + outputdirs[sel_mass_and_age][0]):
                os.makedirs(retrieval_output_dir + outputdirs[sel_mass_and_age][0])

            print('\tCopying %s' % retrieval_output_dir + outputdirs[sel_mass_and_age][0] + '/' + beta_files[sel_beta][0])
            shutil.copyfile(retrieval_input_dir + outputdirs[sel_mass_and_age][0] + '/' + beta_files[sel_beta][0],
                            retrieval_output_dir + outputdirs[sel_mass_and_age][0] + '/' + beta_files[sel_beta][0])

        # if not, interpolate
        else:
            lower_b = b_array[b_array < beta].max()
            upper_b = b_array[b_array > beta].min()

            beta_interpolate(
                input_dir=retrieval_input_dir, 
                output_dir=retrieval_output_dir, 
                current_sub_dir=outputdirs[sel_mass_and_age][0],
                num_interp_points=-99,
                sim_beta_files = [beta_files[b_array == lower_b][0], beta_files[b_array == upper_b][0]],
                beta_arr = np.array([beta])
                )


    # else check if already have this mass
    elif len(sel_mass) != 0:
        print('#2')
        # don't have age b/c checked above
        # get age neighbors
        t_array_at_this_mass = t_array[sel_mass]
        lower_t = t_array_at_this_mass[t_array_at_this_mass < age].max()
        upper_t = t_array_at_this_mass[t_array_at_this_mass > age].min()  

        print('have mass; age neighbors:')
        print(outputdirs[(m_array == mass) & (t_array == lower_t)][0])
        print(outputdirs[(m_array == mass) & (t_array == upper_t)][0])
        
        # interpolate in age
        age_interpolate(
            input_dir=retrieval_input_dir,
            output_dir=retrieval_scratch_dir, 
            mass_string=outputdirs[m_array == mass][0].split('_')[0], 
            t1=outputdirs[(m_array == mass) & (t_array == lower_t)][0].split('_')[1],
            t2=outputdirs[(m_array == mass) & (t_array == upper_t)][0].split('_')[1],
            num_interp_points=-99,
            age_arr=np.array([age])
            )

        # this needs to match how we assign dir name in age_interpolate
        new_dir = outputdirs[m_array == mass][0].split('_')[0] + '_t' + str(round(age, 3))[:5]

        # get betas for this dir
        beta_files = [f for f in os.listdir(retrieval_scratch_dir + new_dir) if not f.startswith('.')]
        b_array = [float(f.split('.dat')[0]) for f in beta_files]

        beta_files = np.array(beta_files)
        b_array = np.array(b_array)

        # check if beta is outside range for given directory
        if not (min(b_array) <= beta <= max(b_array)):
            print('ERROR: beta outside range for this mass and age, which is', [min(b_array),max(b_array)])
            exit()

        # check if already have beta
        sel_beta = np.where(b_array == beta)[0]
        if len(sel_beta) != 0:
            # copy existing beta to retrieval
            if not os.path.exists(retrieval_output_dir + new_dir):
                os.makedirs(retrieval_output_dir + new_dir)

            shutil.copyfile(retrieval_scratch_dir + new_dir + '/' + beta_files[sel_beta][0],
                            retrieval_output_dir + new_dir + '/' + beta_files[sel_beta][0])

        # if not, interpolate
        else:
            lower_b = b_array[b_array < beta].max()
            upper_b = b_array[b_array > beta].min()

            beta_interpolate(
                input_dir=retrieval_scratch_dir, 
                output_dir=retrieval_output_dir, 
                current_sub_dir=new_dir,
                num_interp_points=-99,
                sim_beta_files = [beta_files[b_array == lower_b][0], beta_files[b_array == upper_b][0]],
                beta_arr = np.array([beta])
                )
        

    # else find mass neighbors
    else:
        print('#3')
        lower_m = m_array[m_array < mass].max()
        upper_m = m_array[m_array > mass].min()

        # for lower mass
        t_array_lower_mass = t_array[m_array == lower_m]

        # check if have this age
        if len(t_array_lower_mass[t_array_lower_mass == age]) != 0:
            print('lower mass neighbor:')
            print(outputdirs[(m_array == lower_m) & (t_array == age)][0])

            # TODO this is inefficient
            # copy to retrieval_scratch
            if not os.path.exists(retrieval_scratch_dir):
                os.makedirs(retrieval_scratch_dir)
            # todo doesn't work if already exists
            shutil.copytree(retrieval_input_dir + outputdirs[(m_array == lower_m) & (t_array == age)][0], 
                            retrieval_scratch_dir + outputdirs[(m_array == lower_m) & (t_array == age)][0])

        # if not, find lower mass neighbors
        else:
            lower_m_lower_t = t_array_lower_mass[t_array_lower_mass < age].max()
            lower_m_upper_t = t_array_lower_mass[t_array_lower_mass > age].min()
            print('lower mass neighbors:')
            print(outputdirs[(m_array == lower_m) & (t_array == lower_m_lower_t)][0])
            print(outputdirs[(m_array == lower_m) & (t_array == lower_m_upper_t)][0])
            
            # interpolate in age
            age_interpolate(
                input_dir=retrieval_input_dir,
                output_dir=retrieval_scratch_dir, 
                mass_string=outputdirs[m_array == lower_m][0].split('_')[0], 
                t1=outputdirs[(m_array == lower_m) & (t_array == lower_m_lower_t)][0].split('_')[1],
                t2=outputdirs[(m_array == lower_m) & (t_array == lower_m_upper_t)][0].split('_')[1],
                num_interp_points=-99,
                age_arr=np.array([age])
                )

        # for upper mass
        t_array_upper_mass = t_array[m_array == upper_m]

        # check if have this age
        if len(t_array_upper_mass[t_array_upper_mass == age]) != 0:
            print('upper mass neighbor:')
            print(outputdirs[(m_array == upper_m) & (t_array == age)][0])

            # TODO this is inefficient
            # copy to retrieval_scratch
            if not os.path.exists(retrieval_scratch_dir):
                os.makedirs(retrieval_scratch_dir)
            # todo doesn't work if already exists
            shutil.copytree(retrieval_input_dir + outputdirs[(m_array == upper_m) & (t_array == age)][0], 
                            retrieval_scratch_dir + outputdirs[(m_array == upper_m) & (t_array == age)][0])

        # if not, find lower mass neighbors
        else:
            upper_m_lower_t = t_array_upper_mass[t_array_upper_mass < age].max()
            upper_m_upper_t = t_array_upper_mass[t_array_upper_mass > age].min()
            print('upper mass neighbors:')
            print(outputdirs[(m_array == upper_m) & (t_array == upper_m_lower_t)][0])
            print(outputdirs[(m_array == upper_m) & (t_array == upper_m_upper_t)][0])

            # interpolate in age
            age_interpolate(
                input_dir=retrieval_input_dir,
                output_dir=retrieval_scratch_dir, 
                mass_string=outputdirs[m_array == upper_m][0].split('_')[0], 
                t1=outputdirs[(m_array == upper_m) & (t_array == upper_m_lower_t)][0].split('_')[1],
                t2=outputdirs[(m_array == upper_m) & (t_array == upper_m_upper_t)][0].split('_')[1],
                num_interp_points=-99,
                age_arr=np.array([age])
                )

        # interpolate in mass
        mass_interpolate(
            input_dir=retrieval_scratch_dir, 
            output_dir=retrieval_scratch_dir, 
            # this needs to match how we assign dir name in age_interpolate
            age_string='t'+str(round(age, 3))[:5], 
            m1=outputdirs[m_array == lower_m][0].split('_')[0],
            m2=outputdirs[m_array == upper_m][0].split('_')[0], 
            num_interp_points=-99,
            mass_arr=np.array([mass])
            )

        # this needs to match how we assign dir names in mass_interpolate and age_interpolate
        new_dir = 'm' + str(round(mass, 3))[:5] + '_t' + str(round(age, 3))[:5]

        # get betas for this dir
        beta_files = [f for f in os.listdir(retrieval_scratch_dir + new_dir) if not f.startswith('.')]
        b_array = [float(f.split('.dat')[0]) for f in beta_files]

        beta_files = np.array(beta_files)
        b_array = np.array(b_array)

        # check if beta is outside range for given directory
        if not (min(b_array) <= beta <= max(b_array)):
            print('ERROR: beta outside range for this mass and age, which is', [min(b_array),max(b_array)])
            exit()

        # check if already have beta
        sel_beta = np.where(b_array == beta)[0]
        if len(sel_beta) != 0:
            # copy existing beta to retrieval
            if not os.path.exists(retrieval_output_dir + new_dir):
                os.makedirs(retrieval_output_dir + new_dir)

            shutil.copyfile(retrieval_scratch_dir + new_dir + '/' + beta_files[sel_beta][0],
                            retrieval_output_dir + new_dir + '/' + beta_files[sel_beta][0])

        # if not, interpolate
        else:
            lower_b = b_array[b_array < beta].max()
            upper_b = b_array[b_array > beta].min()

            beta_interpolate(
                input_dir=retrieval_scratch_dir, 
                output_dir=retrieval_output_dir, 
                current_sub_dir=new_dir,
                num_interp_points=-99,
                sim_beta_files = [beta_files[b_array == lower_b][0], beta_files[b_array == upper_b][0]],
                beta_arr = np.array([beta])
                )
            

    # delete scratch directory tree
    # comment out for debug
    if os.path.exists(retrieval_scratch_dir):
        shutil.rmtree(retrieval_scratch_dir)