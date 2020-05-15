import pointLISA 
import os
import matplotlib.pyplot as plt

settings_input = os.getcwd()+'/settings.txt' # Reads settings from this file
orbit_file = '/home/ester/git/synthlisa/orbits/new/Hallion_2pt5G_orbits_pos_uniquedays_timestep_days_scale_1000.txt'

# Obtain STAT object (data, without pointing)
data = pointLISA.static.STAT(settings=settings_input,orbit_file=orbit_file)

# Obtain pointing (aim)
aim = pointLISA.run_din.get_pointing(data) # aim object containing the function for all pointing angles

# Obtain output (out)
i=1 #SC, can be a list of multiple SC
mode='center' # either 'center' (center of aperture), 'mean' (men over ~Nbins**2 of points). 'var' (normalized variance over multiple point) or 'mean_var' (both mean and var
option='both' # either 'function (returns function), 'sampled' (returns sampled values) or both
side='l' #side, can be 'l', 'r', ['l','r'] or 'all'
rets = [['tele_ang','PAAM_ang','offset'],['xoff','yoff','zoff','angx_wf_send','angy_wf_send','angx_wf_rec','angy_wf_rec'],['alpha','Ival']]
t = data.t_all[40]

out_all=[]
for ret in rets:
    out_all.append(pointLISA.calc.values(aim,i,t,side,mode='send',ret=ret))
    for r in ret:
        print(r,getattr(out_all[-1],r))
