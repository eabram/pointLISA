import pointLISA
import os
import matplotlib.pyplot as plt
import numpy as np
import math

settings = os.getcwd()+'/settings/settings_2.txt'
orbit_file = os.getcwd()+'/orbits/Hallion_2pt5G_orbits_pos_uniquedays_timestep_days_scale_1000.txt'

# Obtaining STAT object (static)
data = pointLISA.static.STAT(settings=settings,orbit_file=orbit_file,length_calc=2) # STAT object with length 2 datapoints

# Obtaining AIM object (dynamic)
aim = pointLISA.AIM.AIM(data,tele_control='SS')

# Obtaining figures
t_plot = np.linspace(aim.data.t_all[0],aim.data.t_all[-1],100*len(aim.data.t_all)+1)
SC=1 # Spacecraft number
side='l'
mode='send'
cases = ['tele_send','beam_send','offset_send','xoff','yoff','angx_wf_rec','angy_wf_rec','alpha','Ival']
scales_list = [1e6,1e6,1e6,1,1,1e6,1e6,1e6,1.0]
units_list=['Angle (murad)','Angle (murad)','Angle (murad)','Distance (m)','Distance (m)','Angle (murad)','Angle (murad)','Angle (murad)','Intencity (W/m2)']
scales={}
units={}
for j in range(0,len(cases)):
    scales[cases[j]]=scales_list[j]
    units[cases[j]]=units_list[j]

ret = {}
ret_fast = {}
A = [pointLISA.output.get_output(aim,SC,t,side,mode,cases) for t in t_plot]

for case in cases:
    print(case)
    ret[case] = []
    for a in A:
        ret[case].append(getattr(a,case))

f,ax = plt.subplots(len(ret.keys()),figsize=(21*2,14*len(ret.keys())))
count=0
for k in ret.keys():
    ax[count].plot(t_plot/pointLISA.day2sec,np.array(ret[k])*scales[k])
    ax[count].set_title(k)
    ax[count].set_xlabel('Time (days}')
    ax[count].set_ylabel(units[k])
    count = count+1

f.subplots_adjust(left=0.2, bottom=None, right=None, top=None, wspace=0.3, hspace=2.0)
f.show()





