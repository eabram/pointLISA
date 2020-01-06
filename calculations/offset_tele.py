import pointLISA 
from pointLISA import *
import os

def get_offset(direct=False,length_calc='all'):
    if direct==False:
        direct = os.getcwd()+'/orbit_analysis/'

    if not os.path.exists(direct):
        os.makedirs(direct)

    relativistic=True

    # Options
    set_extr={}
    set_extr['offset_tele'] = 0

    data = pointLISA.run_stat.do_run(length_calc=length_calc,relativistic=False)
    # AIM
    aim0 = pointLISA.run_din.get_pointing(data,tele_control='no_control',PAAM_control='no_control',sampled=False)
    aim = pointLISA.run_din.get_pointing(data,tele_control='full_control',PAAM_control='full_control',aim0=aim0,aim_old=aim0,option_tele='center',option_PAAM='center',sampled=False,set_extr=set_extr)

    #include=['angx_wf_send']
    include=['tele_ang','PAAM_ang']
    side=['l','r']
    i_sel='all'
    mode='center'
    option='sampled'
    out = output.OUTPUT(aim)
    ret = out.make_functions(include=include,option=option,i=i_sel,side=side,mode=mode)

    read_write.write(ret[1],aim.aimset,include=include,direct=direct)
    tele_offset=aim.offset_tele

    for k1 in tele_offset.keys():
        for k2 in tele_offset[k1].keys():
            tele_offset[k1][k2] = -np.nanmean(getattr(getattr(getattr(ret[1],k1),'i'+str(k2)),'angx_wf_send')[0][1])

    offset_file = read_write.write(False,False,direct=direct+'/offset_tele/'+str(length_calc)+'_days/',offset=tele_offset)
    f = open(offset_file,'r')
    set_extr_new = set_extr
    for line in f:
        set_extr_new['offset_tele'] = yaml.load(line)
        break

    print("Obtained telescope offset:")
    print(set_extr_new)

    f.close()
    
    return offset_file,set_extr_new['offset_tele']


