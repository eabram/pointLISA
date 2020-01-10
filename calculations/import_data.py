import pointLISA
from pointLISA import *
import os
import matplotlib.pyplot as plt

def import_data_calc(direct_list):
    f_list=[]
    for direct in direct_list:
        for (dirpath, dirnames, filenames) in os.walk(direct):
            for f in filenames:
                f_list.append(dirpath+'/'+f.split('/')[-1])

    ret=utils.Object
    options=[]
    settings=[]
    for f in f_list:
        ret_new, options_new, settings_new = pointLISA.read_write.read_output(f,ret=ret)
        ret = ret_new
        del ret_new
        options.append(options_new)
        settings.append(settings_new)

    return ret, options, settings

def get_plot(SC,side,ret):
    val_all = {}
    t_all={}

    A = getattr(getattr(ret,side),'i'+str(SC))
    params = A.__dict__.keys()
    for j in params:
        [t,val] = getattr(A,j)[0]
        try:
            val_all[j].append(val)
            t_all[j].append(t)
        except KeyError:
            val_all[j] = val
            t_all[j] = t
    s = 5
    f,ax = plt.subplots(len(params),1,figsize=(5,5*len(params)))
    for j in range(0,len(params)):
        ax[j].plot(np.array(t_all[params[j]])/parameters.day2sec,val_all[params[j]])
        ax[j].set_xlabel('Time (days)')
        ax[j].set_title(params[j])

    return f
