import os
import imports
import utils
import numpy as np

filename = ''
#read_max = 'all'
scale = 'Default'
method = 'fsolve'
new_folder = True
calc_method= 'Abram'
dir_savefig= os.getcwd() +'/' # The directory where the figures will be saved. If False, it will be in the current working directory
noise_check=False
home='/home/ester/git/synthlisa/' # Home directory
directory_imp= False
num_back= 0
dir_orbits= '/home/ester/git/synthlisa/orbits/' # Folder with orbit files
length_calc= 'all' # Length of number of imported datapoints of orbit files. 'all' is also possible
dir_extr= 'zzzAbram_no_abberation' # This will be added to the folder name of the figures
timeunit='Default' # The timeunit of the plots (['minutes'],['days']['years'])
LISA_opt='cache' # If a LISA object from syntheticLISA will be used for further calculations (not sure if it works properly if this False)
arm_influence= True # Set True to consider the travel time of the photons when calculating the nominal armlengths
tstep=False
delay=True #'Not ahead' or False
method='fsolve' # Method used to solve the equation for the photon traveling time
valorfunc='Function' #
select='Hallion' # Select which orbit files will be imported ('all' is all)
aberration=False
delay= True
relativistic= True
test_calc = False
hstep=100 #Time step for calculating the velocity (over hstep seconds average)
putp_mode='sampled' # Or 'LISA'

### AIM settings
aimset = utils.Object()
aimset.inp = False
aimset.tele_control='no_control'
aimset.PAAM_control='no_control'
aimset.tele_ang_extra=False #inplement later
aimset.PAAM_ang_extra=False
aimset.init=False
aimset.sampled=False
aimset.aim_old=False
aimset.aim0=False
aimset.option_tele='center'
aimset.option_PAAM='center'
aimset.optimize_PAAM = 'yoff'
#aimset.optimize_PAAM_value=np.float64(0.0)
aimset.optimize_PAAM_value=np.float64(0.0)
aimset.optimize_PAAM_margin=1000.0
aimset.offset_tele='read'
aimset.tele_method_solve='iter'
aimset.PAAM_method_solve='iter'
aimset.sample_speed =1
aimset.width = 30000.0
aimset.value_center = 0.0
aimset.value_wavefront = 0.0

### Limits/accuracies
aimset.limits = utils.Object()
aimset.limits.xoff = np.float64(1.0e-9) #tele center
aimset.limits.yoff = np.float64(1.0e-9) #PAAM center
aimset.limits.angx = np.float64(1.0e-9) #tele wavefront
aimset.limits.angy = np.float64(1.0e-9) #PAAM wavefront

del os, imports, utils, np
