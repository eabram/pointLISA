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
stat.relativistic= False #Using relativistic calculations (False = classical)
stat.test_calc = False #If True, STAT object would not be written to an object
stat.hstep=0.01 #Time step for calculating the velocity (over hstep seconds average)
stat.putp_mode='pointLISA' # 'interp1d' or 'pointLISA' or 'LISA'

stat.test_COM_effect = False

### AIM settings
aimset = utils.Object()
aimset.inp = False #NN
aimset.tele_control='no_control' # The telescope control method, either 'no_control' (no telescope actuation), 'full_control' (continuously actuating) or 'SS' (Step-and-Stare)
aimset.PAAM_control='no_control' # The one axis PAAM control (inplane) method, either 'no_control' (no telescope actuation), 'full_control' (continuously actuating) or 'SS' (Step-and-Stare). 
aimset.PAAMin_control='no_control' # For a dual axis PAAM the inplane control method, either 'no_control' (no telescope actuation), 'full_control' (continuously actuating) or 'SS' (Step-and-Stare)
aimset.PAAMout_control='no_control' # For a dual axis PAAM the outplane control method, either 'no_control' (no telescope actuation), 'full_control' (continuously actuating) or 'SS' (Step-and-Stare)
aimset.tele_ang_extra=False # If True an offset angle is introduced in the telescoop pointing
aimset.PAAM_ang_extra=False # It True an offset angle is introduces in the PAAM pointing
aimset.offset_tele='0' 
aimset.sampled=False # If True: All coordinates are sampled (and interpolated) and used for further calculations which is faster than when this is set to False, which is slower but more accurate
aimset.aim_old=False 
aimset.aim0=False
aimset.option_tele='center' # The telescope pointing control method either 'center' or wavefront'. The telescope will be pointed to 'center' by aiming the bemline on the center of the opposite aperture and to 'wavefront' by minimizing the angle between the receiveing wavefront and its aperture 
aimset.option_PAAM='center' # The PAAM pointing control method either 'center' or wavefront'. The PAAM will be pointed to 'center' by aiming the bemline on the center of the opposite aperture and to 'wavefront' by minimizing the angle between the receiveing wavefront and its aperture (on the opposite telescope) 
aimset.optimize_PAAM = 'yoff' # The parameter of the PAAM control which is being optimized
aimset.optimize_PAAM_value=np.float64(0.0) # At which aimset.optimize_PAAM has to be optimized for
aimset.optimize_PAAM_margin=10.0 # The optimization margin of aimset.optimize_PAAM
aimset.tele_method_solve='iter' # How, when setting aimset.option_tele to 'wavefront' the optimization is realised (currently only 'iter' (iteration with a certain convergence) is yet implemented
aimset.PAAM_method_solve='iter' # How, when setting aimset.option_PAAM to 'wavefront' the optimization is realised (currently only 'iter' (iteration with a certain convergence) is yet implemented
aimset.sample_speed = 1 # This is either 0 or 1. When in is set to 0 more interpolation points are being sampled then when setting it to 1, however 1 will be faster
aimset.width = 30000.0 # If aimset.tele_method is 'SS' and its method is 'step' (see methods.SS_value()) the maximum of the beamline may be away from the receiving aperture
aimset.power = 1.0e-12 # Minimum received power
aimset.FOV = 8.0e-6 # The Field of View, this value is used when pointing with the Step-and-Stare method
aimset.value_center = 0.0 # If aimset.optimize_tele is set to 'center' this is the value it optimizes for 
aimset.value_wavefront = 0.0 # If aimset.optimize_tele is set to 'wavefront' this is the value it optimizes for 
aimset.PAAM_deg = 1 # Number of rotational axis of PAAM (either 1 or 2)
aimset.tele_SS_scale = 1 # An sclaing parameter for calulating new telescope pointing angles (1 works in can be set to a maximum of 1.9)
aimset.import_file = None # The importet file (with pointing angles), so it does not perform the recalculation

aimset.testSS=False

### Limits/accuracies
aimset.limit_xoff = np.float64(1.0e-9) #tele center
aimset.limit_yoff = np.float64(1.0e-9) #PAAM center
aimset.limit_angx = np.float64(1.0e-9) #tele wavefront
aimset.limit_angy = np.float64(1.0e-9) #PAAM wavefront

del os, imports, utils, np
