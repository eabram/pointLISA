from imports import * 
import copy

### This runner file uses a STAT object (data), calculates the pointig angles and other properties and writes it to a dafile
def get_pointing(data,import_file=None,filename=False,set_din=utils.Object(),aim0=False,aim_old=False,print_on=False,**kwargs):
    from pointLISA import *
    try:
        input_file = data.input_file
    except:
        input_file=None

    aimset0 = utils.get_settings(settings_input=input_file,select='aimset')
    aimset_new={}
    
    for key,value in kwargs.items():
        aimset_new[key] = value

    for k in set_din.__dict__.keys():
        if k not in aimset_new.keys():
            aimset_new[k] = getattr(set_din,k)

    for k in aimset0.__dict__.keys():
        if k not in aimset_new.keys():
            aimset_new[k] = getattr(aimset0,k)
    
    aimset = utils.Object()
    for k in aimset_new.keys():
        if k in aimset0.__dict__.keys():
            setattr(aimset,k,aimset_new[k])
    
    PAAM_deg = aimset.PAAM_deg

    sampled=False
    count=0

    print('check')
    if data.input_file==None:
        if PAAM_deg==1:
            aim = AIM.AIM(import_file=import_file,data=data,setting = aimset,filename=filename,inp=False,aim0=aim0,aim_old=aim0)
            if aim.aimset.tele_control!='SS':
                aim.tele_aim(method=aim.aimset.tele_control,tele_ang_extra=aim.aimset.tele_ang_extra,option=aim.aimset.option_tele)
                option='Default'
            else:
                aim.tele_aim(method='no_control',tele_ang_extra=aim.aimset.tele_ang_extra,option=aim.aimset.option_tele)
                option='SS'
            out = aim.PAAM_aim(method=aim.aimset.PAAM_control,PAAM_ang_extra=aim.aimset.PAAM_ang_extra,option=aim.aimset.option_PAAM)
            if option=='SS':
                aim.tele_aim(method=option,tele_ang_extra=aim.aimset.tele_ang_extra,option=aim.aimset.option_tele)


        
        elif PAAM_deg==2: 
            aim = AIM.AIM(import_file=import_file,data=data,option_tele='center',option_PAAM='center',setting=aimset,init=False,PAAM_deg=2)
            aim.twoPAAM_angles()
            print('Under construction')
            


    else:
        print('tele:')
        print(aimset.tele_control)
        print('PAAM:')
        print(aimset.PAAM_control)

        aim = AIM.AIM(data=data,setting=aimset,filename=filename,inp=False)
        aim.tele_aim(method=aim.aimset.tele_control,tele_ang_extra=aim.aimset.tele_ang_extra,option=aim.aimset.option_tele)
        out = aim.PAAM_aim(method=aim.aimset.PAAM_control,PAAM_ang_extra=aim.aimset.PAAM_ang_extra,option=aim.aimset.option_PAAM)

    return aim
