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

for outputdir in outputdirs:
    m = outputdir.split('_')[0][1:]
    t = outputdir.split('_')[1][1:]

    if m == mass:
        print('already have this mass')
        if t == age:
            print('already have this age')

            betadirs = os.listdir('../output/'+outputdir)
            for betadir in betadirs:
                b = betadir.split('.dat')[0]
                if b == beta:
                    print('already have this beta')
                    
                    print('Your requested dM/dt is here:')
                    print('../output/' + outputdir + '/' + betadir)
