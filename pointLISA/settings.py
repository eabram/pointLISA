import os
import imports
import utils
import numpy as np

# Static settings
stat = utils.Object()
stat.filename = ''
stat.read_max = 'all'
stat.scale = 'Default'
stat.method = 'fsolve'
stat.new_folder = True
stat.calc_method= 'Abram'
stat.dir_savefig= os.getcwd() +'/' # The directory where the figures will be saved. If False, it will be in the current working directory
stat.noise_check=False
stat.home='/home/ester/git/synthlisa/' # Home directory
stat.directory_imp= False
stat.num_back= 0
stat.dir_orbits= '/home/ester/git/synthlisa/orbits/' # Folder with orbit files
stat.length_calc= 'all' # Length of number of imported datapoints of orbit files. 'all' is also possible
stat.dir_extr= 'zzzAbram_no_abberation' # This will be added to the folder name of the figures
stat.timeunit='Default' # The timeunit of the plots (['minutes'],['days']['years'])
stat.LISA_opt='cache' # If a LISA object from syntheticLISA will be used for further calculations (not sure if it works properly if this False)
stat.arm_influence= True # Set True to consider the travel time of the photons when calculating the nominal armlengths
stat.tstep=False
stat.delay=True #'Not ahead' or False
stat.method='fsolve' # Method used to solve the equation for the photon traveling time
stat.valorfunc='Function' #
stat.select='Hallion' # Select which orbit files will be imported ('all' is all)
stat.aberration=False
stat.delay= True
stat.relativistic= True
stat.test_calc = False
stat.hstep=100 #Time step for calculating the velocity (over hstep seconds average)
stat.putp_mode='sampled' # Or 'LISA'

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
#aimset.limits = utils.Object()
aimset.limit_xoff = np.float64(1.0e-9) #tele center
aimset.limit_yoff = np.float64(1.0e-9) #PAAM center
aimset.limit_angx = np.float64(1.0e-9) #tele wavefront
aimset.limit_angy = np.float64(1.0e-9) #PAAM wavefront

del os, imports, utils, np
