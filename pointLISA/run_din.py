from imports import *

def get_pointing(data,filename=False,**kwargs):
    aimset = settings.aimset
    for key,value in kwargs.items():
        setattr(aimset,key,value)

    aim0 = kwargs.pop('aim0',False)
    aim_old = kwargs.pop('aim_old',False)

    sampled=False
    count=0

    print('check')

#    if aimset.PAAM_ang_extra==True:
#        aimset.PAAM_ang_extra = methods.get_extra_ang_mean(aimsetself,'PAAM')

    if aim_old==False:
        tele_l_ang=False
        tele_r_ang=False
        PAAM_l_ang=False
        PAAM_r_ang=False
    else:
        tele_l_ang=aim_old.tele_l_ang
        tele_r_ang=aim_old.tele_r_ang
        PAAM_l_ang=aim_old.beam_l_ang
        PAAM_r_ang=aim_old.beam_r_ang
    angles_old = [tele_l_ang,PAAM_l_ang,tele_r_ang,PAAM_r_ang]

    if aim0==False:
        tele_l_ang0=False
        tele_r_ang0=False
        PAAM_l_ang0=False
        PAAM_r_ang0=False
    else:
        tele_l_ang0=aim0.tele_l_ang
        tele_r_ang0=aim0.tele_r_ang
        PAAM_l_ang0=aim0.beam_l_ang
        PAAM_r_ang0=aim0.beam_r_ang
    angles0 = [tele_l_ang0,PAAM_l_ang0,tele_r_ang0,PAAM_r_ang0]

    aim = AIM.AIM(data=data,init=aimset.init,sampled=aimset.sampled,angles0=angles0,angles_old=angles_old,offset_tele=aimset.offset_tele,settings=aimset,filename=filename,inp=aimset.inp)
    aim.tele_aim(method=aim.aimset.tele_control,tele_ang_extra=aim.aimset.tele_ang_extra,option=aim.aimset.option_tele)
    out = aim.PAAM_aim(method=aim.aimset.PAAM_control,PAAM_ang_extra=aim.aimset.PAAM_ang_extra,option=aim.aimset.option_PAAM)
    
    return out



