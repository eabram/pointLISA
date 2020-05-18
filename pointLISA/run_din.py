from pointLISA import * 

### This runner file uses a STAT object (data), calculates the pointig angles and other properties and writes it to a dafile
def get_pointing(data,import_file=None,filename=False,**kwargs):
    try:
        setting_file = data.settings
    except:
        input_file=None

    aimset = const.get_settings(settings_input=setting_file,select='aimset',kwargs=kwargs)
    
    PAAM_deg = aimset.PAAM_deg

    sampled=False
    count=0

    print('check')
    if data.stat.filename==None:
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
