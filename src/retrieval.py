"""
Aim:
One should be able to specify a (list of) parameters (mass, age, beta) to retrieve a dM/dt curve for.
I think we want to implement this as a quick 3D interpolation for each row in the list.
"""

inputs = [
#mass,     age,    beta
[1,        0.5,    2.2]
]

#for row in inputs:

# todo just try for one first
row = inputs[0]
mass = row[0]
age = row[1]
beta = row[2]

# look if mass is already in output
outputdirs = os.listdir('../output/')

found_mass = False
found_age = False
found_beta = False

for outputdir in outputdirs:
    m = outputdir.split('_')[0][1:]
    t = outputdir.split('_')[1][1:]

    if m == mass:
        print('already have this mass')
        found_mass = True

        if t == age:
            print('already have this age')
            found_age = True

            betadirs = os.listdir('../output/'+outputdir)
            for betadir in betadirs:
                b = betadir.split('.dat')[0]
                if b == beta:
                    print('already have this beta')
                    found_beta = True

                    print('Your requested dM/dt is here:')
                    print('../output/' + outputdir + '/' + betadir)

                    break

if (found_mass == True) & (found_age == True) & (found_beta == True):
    print('We are done.')
    exit()

if (found_mass == True) & (found_age == True) & (found_beta == False):
    print('Found mass, age, not beta.')
    ## interpolate in beta between 2 closest betas


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