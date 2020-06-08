import pointLISA
import os
import matplotlib.pyplot as plt
import numpy as np
import math

settings = os.getcwd()+'/settings/settings_c.txt'
orbit_file = os.getcwd()+'/orbits/Hallion_2pt5G_orbits_pos_uniquedays_timestep_days_scale_1000.txt'

# Obtaining STAT object (static)
data = pointLISA.static.STAT(settings=settings,orbit_file=orbit_file) # STAT object with default values
data_pl = pointLISA.static.STAT(settings=settings,orbit_file=orbit_file,LISA_opt='pointLISA') #STAT object with changed setting (LISA_opt is the interpolation method)

# Obtaining AIM object (dynamic)
aim = pointLISA.AIM.AIM(data)
aim_pl = pointLISA.AIM.AIM(data_pl)

# Obtaining figures
t_plot = aim.data.t_all[5:-10]
SC=1 # Spacecraft number

L = np.array([aim.data.L_sl(SC,t) for t in t_plot]) # Photon traveling time when exeting de left telescope on spacecraft SC
L_pl = np.array([aim_pl.data.L_sl(SC,t) for t in t_plot])

f,ax = plt.subplots(2,1,figsize=(5,10))
ax[0].plot(t_plot/pointLISA.day2sec,L,label='syntheticLISA fit')
ax[0].plot(t_plot/pointLISA.day2sec,L_pl,label='pointLISA fit')
ax[0].set_title('Photon traveling time')
ax[0].set_xlabel('Time (days)')
ax[0].set_ylabel('Time (seconds)')
ax[1].plot(t_plot/pointLISA.day2sec,L-L_pl)
ax[1].set_title('Photon traveling time difference')
ax[1].set_xlabel('Time (days)')
ax[1].set_ylabel('Time (seconds)')
f.subplots_adjust(left=0.2, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3)
f.show()

t_plot = aim.data.t_all[5:20]
side='l'
mode='send'
cases = ['tele_send','beam_send','xoff','yoff','angx_wf_rec','angy_wf_rec']
scales_list = [1e6,1e6,1,1,1e6,1e6]
units_list=['Angle (murad)','Angle (murad)','Distance (m)','Distance (m)','Angle (murad)','Angle (murad)']
scales={}
units={}
for j in range(0,len(cases)):
    scales[cases[j]]=scales_list[j]
    units[cases[j]]=units_list[j]

ret = {}
ret_fast = {}
A = [pointLISA.output.get_output(aim,SC,t,side,mode,cases) for t in t_plot]
B = [pointLISA.output.get_output(aim_pl,SC,t,side,mode,cases) for t in t_plot]

for case in cases:
    print(case)
    ret[case] = []
    ret_fast[case] = []
    for a in A:
        ret[case].append(getattr(a,case))
    for b in B:
        ret_fast[case].append(getattr(b,case))

f,ax = plt.subplots(len(ret.keys()),2,figsize=(21*2,14*len(ret.keys())))
count=0
for k in ret.keys():
    ax[count,0].plot(t_plot/pointLISA.day2sec,np.array(ret[k])*scales[k])
    ax[count,0].plot(t_plot/pointLISA.day2sec,np.array(ret_fast[k])*scales[k])  
    ax[count,0].set_title(k)
    ax[count,0].set_xlabel('Time (days}')
    ax[count,0].set_ylabel(units[k])
    ax[count,1].plot(t_plot/pointLISA.day2sec,(np.array(ret[k])-np.array(ret_fast[k]))*scales[k])  
    ax[count,1].set_title(k)
    ax[count,1].set_xlabel('Time (days}')
    ax[count,1].set_ylabel(units[k])
    
    count = count+1

for j in range(0,len(ax)):
    for k in range(0,len(ax[j])):
        ax[j,k].set_xlabel('Time (days')

f.subplots_adjust(left=0.2, bottom=None, right=None, top=None, wspace=0.3, hspace=1.5)
f.show()





