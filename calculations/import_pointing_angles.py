from pointLISA import *
import os

def get_all(input_file=None,set_stat=utils.Object(),set_din=utils.Object()):
    data = pointLISA.run_stat.do_run(input_file=input_file,set_stat=set_stat)
    aim0= pointLISA.run_din.get_pointing(data,tele_control='no_control',PAAM_control='no_control',sampled=False,tele_ang_extra=False,PAAM_ang_extra=False,set_din=set_din,inp=False,init=True)
    aim= pointLISA.run_din.get_pointing(data,aim0=aim0,aim_old=aim0,sampled=False,set_din=set_din,init=False)

    return aim

