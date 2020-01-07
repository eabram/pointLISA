import os  
import imports
import utils
import numpy as np

# Static settings
stat = utils.Object()
stat.scale = 'Default' #The scaling of the imported coordinates ('Default'=scale is set by the filename or equal to 1)
stat.calc_method= 'Abram' # Use either the 'Abram' method or 'Waluschka' method
stat.home='/home/ester/git/synthlisa/' # Home directory
stat.directory_imp= False # If not False adding a subfolder to sta.home 
stat.dir_orbits= '/home/ester/git/synthlisa/orbits/' # Folder with orbit files
stat.length_calc= 'all' # Length of number of imported datapoints of orbit files. 'all' is also possible
#stat.dir_extr= 'zzzAbram_no_abberation' #NN # This will be added to the folder name of the figures
stat.timeunit='Default' # The timeunit of the plots (['minutes'],['days']['years'])
stat.LISA_opt='cache' # If a LISA object from syntheticLISA will be used for further calculations (not sure if it works properly if this False)
stat.delay=True #'Not ahead' or False
stat.valorfunc='Function' #
stat.select='Hallion' # Select which orbit files will be imported ('all' is all)
stat.aberration=True #Consider the aberration angle of the incoming light 
stat.relativistic= True #Using relativistic calculations (False = classical)
stat.test_calc = False #If True, STAT object would not be written to an object
stat.hstep=100 #Time step for calculating the velocity (over hstep seconds average)
stat.putp_mode='sampled' # 'samped' or 'LISA'

### AIM settings
aimset = utils.Object()
aimset.inp = False
aimset.tele_control='no_control'
aimset.PAAM_control='no_control'
aimset.PAAMin_control='no_control'
aimset.PAAMout_control='no_control'
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
aimset.FOV = 8.0e-6
aimset.value_center = 0.0
aimset.value_wavefront = 0.0
aimset.PAAM_deg = 1 # Number of rotational axis of PAAM (either 1 or 2)
aimset.tele_SS_scale = 1
### Limits/accuracies
#aimset.limits = utils.Object()
aimset.limit_xoff = np.float64(1.0e-9) #tele center
aimset.limit_yoff = np.float64(1.0e-9) #PAAM center
aimset.limit_angx = np.float64(1.0e-9) #tele wavefront
aimset.limit_angy = np.float64(1.0e-9) #PAAM wavefront

del os, imports, utils, np
