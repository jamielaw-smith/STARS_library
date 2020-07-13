import os
from configparser import RawConfigParser
import time as tm
from src.interpolation import beta_interpolate, mass_interpolate, age_interpolate
from src.retrieval import retrieval
import numpy as np

config = RawConfigParser()
config.read(os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'STARS.config')))

class STARS_library:

    def add_options(self, parser=None, usage=None, config=None):
        import optparse

        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option("-r", "--retrieve", nargs=3, type="float",
                          help="Tuple to retrieve single model at command line: {mass} {age} {beta}.")

        parser.add_option("-g", "--retrieve_grid", action="store_true", dest="retrieve_grid", default=False, 
                          help="Flag to retrieve grid of models from file, which is RETRIEVE.par by default.\
                           DEFAULT=False")

        return (parser)

    def initialize(self):

        ## BETA INTERPOLATION ##
        all_input_subdirs = [name for name in os.listdir(self.input_dir)
                             if os.path.isdir(os.path.join(self.input_dir, name))]

        for subdir in all_input_subdirs:
            beta_interpolate(self.input_dir, self.output_dir, subdir, self.num_interp_beta)


        ## MASS INTERPOLATION ##
        current_output_subdirs = [name for name in os.listdir(self.output_dir)
                             if os.path.isdir(os.path.join(self.output_dir, name))]

        zams_key = "t0.0"
        tams_key = "t1.0"
        model_directories_by_time = {
            zams_key: None,
            tams_key: None
        }

        zams_masses = []
        tams_masses = []
        for cod in current_output_subdirs:
            dir_name = cod.split("/")[-1]
            age_key = dir_name.split("_")[-1]
            mass_key = dir_name.split("_")[0]

            # We're currently not handling the t0.57 case in the mass interpolation, so skip
            if age_key == "t0.57":
                continue

            if zams_key in dir_name:
                zams_masses.append(mass_key)
            else:
                tams_masses.append(mass_key)

        model_directories_by_time[zams_key] = sorted(zams_masses)
        model_directories_by_time[tams_key] = sorted(tams_masses)

        for age_string, mass_steps in model_directories_by_time.items():
            for m1, m2 in zip(mass_steps[:-1], mass_steps[1:]):
                mass_interpolate(self.output_dir, age_string, m1, m2, self.num_interp_mass)


        ## AGE INTERPOLATION ##
        model_directories_by_mass = {}
        current_output_subdirs = [name for name in os.listdir(self.output_dir)
                                  if os.path.isdir(os.path.join(self.output_dir, name))]
        for sub_dir in current_output_subdirs:
            mass_str = sub_dir.split("_")[0]

            if mass_str not in model_directories_by_mass:
                if mass_str != "m1.0":
                    model_directories_by_mass[mass_str] = ["t0.0", "t1.0"]
                else:
                    model_directories_by_mass[mass_str] = ["t0.0", "t0.57", "t1.0"]

        for mass_string, age_steps in model_directories_by_mass.items():
            for t1, t2 in zip(age_steps[:-1], age_steps[1:]):
                age_interpolate(self.output_dir, mass_string, t1, t2, self.num_interp_age)

    def __init__(self):

        self.input_dir = config.get('general_settings', 'input_dir')
        self.output_dir = config.get('general_settings', 'output_dir')
        self.retrieval_grid_file = config.get('general_settings', 'retrieval_grid_file')
        self.retrieval_input_dir = config.get('general_settings', 'retrieval_input_dir')
        self.retrieval_scratch_dir = config.get('general_settings', 'retrieval_scratch_dir')
        self.retrieval_output_dir = config.get('general_settings', 'retrieval_output_dir')
        self.retrieval_grid_file = config.get('general_settings', 'retrieval_grid_file')
        self.num_interp_beta = int(config.get('interpolation', 'NUM_BETA_INTERP_POINTS'))
        self.num_interp_mass = int(config.get('interpolation', 'NUM_MASS_INTERP_POINTS'))
        self.num_interp_age = int(config.get('interpolation', 'NUM_AGE_INTERP_POINTS'))

        # Check if initialized, if not, initialize first
        # JLS added this as quick fix b/c below fails if ../output/ doesn't exist
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        
        interpolated_subdirs = [name for name in os.listdir(self.output_dir)
                                if os.path.isdir(os.path.join(self.output_dir, name))]

        if len(interpolated_subdirs) == 0:
            print("Initializing library...")

            t1 = tm.time()
            self.initialize()
            t2 = tm.time()

            print("... Initialize complete. [%0.2f sec]" % (t2 - t1))

    def retrieve(self, mass, age, beta):
        retrieval(mass, age, beta, self.retrieval_input_dir, self.retrieval_scratch_dir, self.retrieval_output_dir)


if __name__ == "__main__":
    useagestring = """python STARS_library.py [options]

    ### To run a single model: ###
    `python STARS_library.py -r {mass} {age} {beta}`
    
    retrieve params:
        mass [M_sun]:                               0.3 - 3.0
        age [fractional; 0 == ZAMS, 1.0 == TAMS]:   0.0 - 1.0
        beta [impact parameter]:                    varies
    
    
    
    ### To run a grid of models: ###
    1) Edit the ./RETRIEVE.par file with parameters you wish:
    
    ex:
    
    mass [M_sun]    age [fractional; 0=ZAMS, 1==TAMS]   beta [r_p/r_t]
    1               0                                   1
    2               0.5                                 3
    ...
    
    2) Run the command:
    `python STARS_library.py -g`
    """
    #start = tm.time()

    stars_lib = STARS_library()
    parser = stars_lib.add_options(usage=useagestring)
    options, args = parser.parse_args()
    stars_lib.options = options


    if options.retrieve and options.retrieve_grid:
        raise Exception("Use only -r OR -g.")

    if options.retrieve:
        mass = options.retrieve[0]
        age = options.retrieve[1]
        beta = options.retrieve[2]
        stars_lib.retrieve(mass, age, beta)

    if options.retrieve_grid:
        mass, age, beta = np.loadtxt(stars_lib.retrieval_grid_file, skiprows=1, unpack=True)
        for m, a, b in zip(mass, age, beta):
            stars_lib.retrieve(m, a, b)

    #end = tm.time()
    #duration = (end - start)
    #print("\n********* start DEBUG ***********")
    #print("STARS_library `retrieve` execution time: %s" % duration)
    #print("********* end DEBUG ***********\n")
