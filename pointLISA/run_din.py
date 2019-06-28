from imports import *
import copy

def get_pointing(data,import_file=None,filename=False,set_din=utils.Object(),aim0=False,aim_old=False,print_on=False,**kwargs):
    from pointLISA import *
    try:
        input_file = data.input_file
    except:
        input_file=None

    aimset_din = settings.aimset
    aimset_new = utils.Object()
    if input_file!=None:
        aimset = read_write.read_options(input_file).__dict__
        for k in aimset_din.__dict__.keys():
            try:
                setattr(aimset_new,k,aimset[k])
            except KeyError:
                setattr(aimset_new,k,aimset_din.__dict__[k])
                pass

        del aimset, aimset_din
        aimset = copy.copy(aimset_new)
        del aimset_new
    else:
        aimset = copy.copy(aimset_din)
            #print(dirpath.split(source_folder)[-1]+'/')
            #print(dirpath.split(source_folder)[-1]+'/')
        del aimset_din

#    for k in aimset.__dict__.keys():
#        print(k,getattr(aimset,k))

    if type(set_din)==dict:
        new = utils.Object
        for k in set_din.keys():
            if '__' not in k:
                setattr(new,k,set_din[k])
        del set_din
        set_din = copy.copy(new)
        del new

    for k in set_din.__dict__.keys():
        if k in aimset.__dict__.keys():
            setattr(aimset,k,getattr(set_din,k))
            if print_on:
                print('Adjust setting: '+k+' = '+str(getattr(aimset,k)))
        else:
            print(str(k)+' is not a used option')
    
    for key,value in kwargs.items():
        setattr(aimset,key,value)
        setattr(set_din,key,value)
        if print_on:
            print('Adjust setting: '+key+' = '+str(getattr(aimset,key)))


    sampled=False
    count=0

    print('check')

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

    

    if data.input_file==None:    
        aim = AIM.AIM(data=data,setting = set_din,init=aimset.init,sampled=aimset.sampled,angles0=angles0,angles_old=angles_old,offset_tele=aimset.offset_tele,settings=aimset,filename=filename,inp=False)
        aim.tele_aim(method=aim.aimset.tele_control,tele_ang_extra=aim.aimset.tele_ang_extra,option=aim.aimset.option_tele)
        out = aim.PAAM_aim(method=aim.aimset.PAAM_control,PAAM_ang_extra=aim.aimset.PAAM_ang_extra,option=aim.aimset.option_PAAM)
    else:
        print('tele:')
        print(aimset.tele_control)
        print('PAAM:')
        print(aimset.PAAM_control)

        aim = AIM.AIM(data=data,init=aimset.init,sampled=aimset.sampled,angles0=angles0,angles_old=angles_old,offset_tele=aimset.offset_tele,settings=aimset,filename=filename,inp=False,aimset_read=aimset)
        aim.tele_aim(method=aimset.tele_control,tele_ang_extra=aimset.tele_ang_extra,option=aimset.option_tele)

        out = aim.PAAM_aim(method=aimset.PAAM_control,PAAM_ang_extra=aimset.PAAM_ang_extra,option=aimset.option_PAAM)

    
    #return out

    return aim



