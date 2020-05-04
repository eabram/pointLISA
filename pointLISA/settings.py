from pointLISA import *

# Static settings
#stat = utils.Object()
stat_scale = 'Default' #The scaling of the imported coordinates ('Default'=scale is set by the filename or equal to 1)
stat_calc_method= 'Abram' # Use either the 'Abram' method or 'Waluschka' method
stat_home='/home/ester/git/synthlisa/' # Home directory
stat_directory_imp= False # If not False adding a subfolder to sta.home 
stat_dir_orbits= '/home/ester/git/synthlisa/orbits/' # Folder with orbit files
stat_length_calc= 400 # Length of number of imported datapoints of orbit files. 'all' is also possible
#sta_.dir_extr= 'zzzAbram_no_abberation' #NN # This will be added to the folder name of the figures
stat_timeunit='Default' # The timeunit of the plots (['minutes'],['days']['years'])
stat_LISA_opt='cache' # If a LISA object from syntheticLISA will be used for further calculations (not sure if it works properly if this False)
stat_delay=True #'Not ahead' or False
stat_valorfunc='Function' #
stat_select='Hallion' # Select which orbit files will be imported ('all' is all)
stat_aberration=True #Consider the aberration angle of the incoming light 
stat_relativistic= False #Using relativistic calculations (False = classical)
stat_test_calc = False #If True, STAT object would not be written to an object
stat_hstep=0.01 #Time step for calculating the velocity (over hstep seconds average)
stat_putp_mode='pointLISA' # 'interp1d' or 'pointLISA' or 'LISA'

stat_tidal=0 # Include tidal effects in orbital function
stat_orbit_function=False # Obtain orbital function instead of importing one
stat_test_COM_effect = False

### AIM settings
#aimset = utils.Object()
aimset_inp = False #NN
aimset_tele_control='no_control' # The telescope control method, either 'no_control' (no telescope actuation), 'full_control' (continuously actuating) or 'SS' (Step-and-Stare)
aimset_PAAM_control='no_control' # The one axis PAAM control (inplane) method, either 'no_control' (no telescope actuation), 'full_control' (continuously actuating) or 'SS' (Step-and-Stare). 
aimset_PAAMin_control='no_control' # For a dual axis PAAM the inplane control method, either 'no_control' (no telescope actuation), 'full_control' (continuously actuating) or 'SS' (Step-and-Stare)
aimset_PAAMout_control='no_control' # For a dual axis PAAM the outplane control method, either 'no_control' (no telescope actuation), 'full_control' (continuously actuating) or 'SS' (Step-and-Stare)
aimset_tele_ang_extra=False # If True an offset angle is introduced in the telescoop pointing
aimset_PAAM_ang_extra=False # It True an offset angle is introduces in the PAAM pointing
aimset_offset_tele='0' 
aimset_sampled=False # If True: All coordinates are sampled (and interpolated) and used for further calculations which is faster than when this is set to False, which is slower but more accurate
aimset_aim_old=False 
aimset_aim0=False
aimset_option_tele='center' # The telescope pointing control method either 'center' or wavefront'. The telescope will be pointed to 'center' by aiming the bemline on the center of the opposite aperture and to 'wavefront' by minimizing the angle between the receiveing wavefront and its aperture 
aimset_option_PAAM='center' # The PAAM pointing control method either 'center' or wavefront'. The PAAM will be pointed to 'center' by aiming the bemline on the center of the opposite aperture and to 'wavefront' by minimizing the angle between the receiveing wavefront and its aperture (on the opposite telescope) 
aimset_optimize_PAAM = 'yoff' # The parameter of the PAAM control which is being optimized
aimset_optimize_PAAM_value=np.float64(0.0) # At which aimset.optimize_PAAM has to be optimized for
aimset_optimize_PAAM_margin=10.0 # The optimization margin of aimset.optimize_PAAM
aimset_tele_method_solve='iter' # How, when setting aimset.option_tele to 'wavefront' the optimization is realised (currently only 'iter' (iteration with a certain convergence) is yet implemented
aimset_PAAM_method_solve='iter' # How, when setting aimset.option_PAAM to 'wavefront' the optimization is realised (currently only 'iter' (iteration with a certain convergence) is yet implemented
aimset_sample_speed = 1 # This is either 0 or 1. When in is set to 0 more interpolation points are being sampled then when setting it to 1, however 1 will be faster
aimset_width = 30000.0 # If aimset.tele_method is 'SS' and its method is 'step' (see methods.SS_value()) the maximum of the beamline may be away from the receiving aperture
aimset_power = 1.0e-12 # Minimum received power
aimset_FOV = 8.0e-6 # The Field of View, this value is used when pointing with the Step-and-Stare method
aimset_value_center = 0.0 # If aimset.optimize_tele is set to 'center' this is the value it optimizes for 
aimset_value_wavefront = 0.0 # If aimset.optimize_tele is set to 'wavefront' this is the value it optimizes for 
aimset_PAAM_deg = 1 # Number of rotational axis of PAAM (either 1 or 2)
aimset_tele_SS_scale = 1 # An sclaing parameter for calulating new telescope pointing angles (1 works in can be set to a maximum of 1.9)
aimset_import_file = None # The importet file (with pointing angles), so it does not perform the recalculation

aimset_testSS=False

### Limits/accuracies
aimset_limit_xoff = np.float64(1.0e-9) #tele center
aimset_limit_yoff = np.float64(1.0e-9) #PAAM center
aimset_limit_angx = np.float64(1.0e-9) #tele wavefront
aimset_limit_angy = np.float64(1.0e-9) #PAAM wavefront
