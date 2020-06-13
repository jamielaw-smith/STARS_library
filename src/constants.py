import numpy as np
import astropy.constants as c

C_CGS = c.c.cgs.value
FOUR_PI = 4*np.pi
M_SUN_CGS = c.M_sun.cgs.value
G = c.G.cgs.value  # 6.67259e-8 cm3 g-1 s-2

model_dir_formatter = "{}_{}"