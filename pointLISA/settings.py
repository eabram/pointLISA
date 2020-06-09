from pointLISA import *

# Static settings
stat_scale = 'Default' #The scaling of the imported coordinates ('Default'=scale is set by the filename or equal to 1)
stat_length_calc= 400 # Length of number of imported datapoints of orbit files. 'all' is also possible
stat_timeunit='Default' # The timeunit of the plots (['minutes'],['days']['years'])
stat_interpolation_method='cache' # If a LISA object from syntheticLISA will be used for further calculations (not sure if it works properly if this False) 'interp1d' or 'pointLISA' or 'LISA'
stat_aberration=True #Consider the aberration angle of the incoming light 
stat_relativistic= False #Using relativistic calculations (False = classical)
stat_hstep=0.01 #Time step for calculating the velocity (over hstep seconds average)

stat_tidal=0 # Include tidal effects in orbital function
stat_test_COM_effect = False

### AIM settings
allignmentset_tele_control='no_control' # The telescope control method, either 'no_control' (no telescope actuation), 'full_control' (continuously actuating) or 'SS' (Step-and-Stare)
allignmentset_PAAM_control='no_control' # The one axis PAAM control (inplane) method, either 'no_control' (no telescope actuation), 'full_control' (continuously actuating) or 'SS' (Step-and-Stare). 
allignmentset_PAAMin_control='no_control' # For a dual axis PAAM the inplane control method, either 'no_control' (no telescope actuation), 'full_control' (continuously actuating) or 'SS' (Step-and-Stare)
allignmentset_PAAMout_control='no_control' # For a dual axis PAAM the outplane control method, either 'no_control' (no telescope actuation), 'full_control' (continuously actuating) or 'SS' (Step-and-Stare)
allignmentset_offset_tele=True 
allignmentset_option_tele='center' # The telescope pointing control method either 'center' or wavefront'. The telescope will be pointed to 'center' by aiming the bemline on the center of the opposite aperture and to 'wavefront' by minimizing the angle between the receiveing wavefront and its aperture 
allignmentset_option_PAAM='center' # The PAAM pointing control method either 'center' or wavefront'. The PAAM will be pointed to 'center' by aiming the bemline on the center of the opposite aperture and to 'wavefront' by minimizing the angle between the receiveing wavefront and its aperture (on the opposite telescope) 
allignmentset_PAAM_deg = 1 # Number of rotational axis of PAAM (either 1 or 2)
allignmentset_import_file = None # The imported file (with pointing angles), so it does not perform the recalculation
allignmentset_alignment_object = None
allignmentset_solve_method='solve'
