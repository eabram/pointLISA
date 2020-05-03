from pointLISA import * 
import os

def get_SS(PAAM_control='full_control',direct=False,set_stat=utils.Object(),set_din = utils.Object()):
    if direct==False:
        direct = os.getcwd()+'/orbit_analysis/'

    if not os.path.exists(direct):
        os.makedirs(direct)

    length_calc=20
    relativistic=True

    # Options
    set_stat.length_calc=20
    set_stat.relativistic=True
    set_din.offset_tele='read'

    data = run_stat.do_run(set_stat=set_stat)

    # AIM
    aim0 = run_din.get_pointing(data,tele_control='no_control',PAAM_control='no_control',sampled=False,set_sin=set_din)
    aim = run_din.get_pointing(data,tele_control='SS',PAAM_control=PAAM_control,aim0=aim0,aim_old=aim0,sampled=True,set_dinextr=set_din)






    include=['tele_l_ang','tele_r_ang','t_l_adjust','t_r_adjust','beam_l_ang','beam_r_ang']
    side=['l','r']
    i_sel='all'
    mode='center'
    option='sampled'
    out = output.OUTPUT(aim)
    ret = out.make_functions(include=include,option=option,i=i_sel,side=side,mode=mode)

    read_write.write(ret[1],aim.aimset,include=include,direct=direct)

    print("Obtained SS")

    return ret

