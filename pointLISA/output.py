from pointLISA import *

def get_output(aim,i_all,t_all,side_all,mode,cases,state='cases',return_format='normal'): # if return_format is 'write'
    options={}
    for k in locals().keys():
        options[k] = locals()[k]

    def get_output_time(aim,i,t,side,mode,cases,state='cases'):
        pos = utils.Object()
        attr = aim.get_value(i,t,side,mode=mode,ret='all')
        for k in attr.keys():
            setattr(pos,k,attr[k])
        setattr(pos,'aim',aim)

        if type(cases)!=list:
            cases=[cases]
        for case in cases:
            try:
                pos = calc_value(pos,case)
            except UnboundLocalError,e:
                print('Could not obtain '+case)
                setattr(pos,case,None)
        if state=='all':
            return pos
        elif state=='cases':
            pos_ret = utils.Object()
            for case in cases:
                setattr(pos_ret,case,getattr(pos,case))
            return pos_ret
        else:
            raise ValueError('Please select state= "all" or "cases"')
    
    try:
        len(side_all)
    except TypeError:
        side_all=[side_all]
    try:
        len(i_all)
    except TypeError:
        i_all = [i_all]

    ret_full = utils.Object()
    for side in side_all:
        setattr(ret_full,side,utils.Object)
        for i in i_all:
            try:
                len(t_all)
                if len(t_all)==1:
                    raise TypeError()
                results=[]
                ret = {}
                for t in t_all:
                    results.append(get_output_time(aim,i,t,side,mode,cases,state='cases'))
                for case in cases:
                    ret[case] = []
                    print(case)
                    for j in range(0,len(t_all)):
                        ret[case].append(getattr(results[j],case))
                    ret[case] = np.array(ret[case])
            except TypeError:
                ret = get_output_time(aim,i,t_all,side,mode,cases,state='cases').__dict__
            ret['t']=np.array(t_all)
            
            setattr(getattr(ret_full,side),'SC'+str(i),ret)
    
    if return_format=='normal' and len(side_all)==1 and len(i_all)==1:
        ret_all = utils.Object
        setattr(ret_all,'options',options)
        setattr(ret_all,'ret',ret)
    elif return_format=='write' or (len(side_all)!=1 and len(i_1ll)!=1):
        ret_all = ret_full
        setattr(ret_all,'options',options)

    return ret_all

def calc_value(pos,case):
    done=False
    case_all=[case]
    loop=True
    while loop is True:
        [pos,done,case_new] = loop_over_cases(pos,case_all[-1],done=done)
        if case_new==case_all[0] and done==True:
            loop = False
            break  
        if done==True:
            case_all.remove(case_new)
        else:
            if case_new not in case_all:
                case_all.append(case_new)

    return pos

def loop_over_cases(pos,case,done=False):
    try:
        pos = get_case(pos,case)
        done = True
        #print('Obtained '+case,getattr(pos,case))
    except AttributeError,e:
        case = str(e).split('has no attribute ')[-1][1:-1]
        done = False
        #print(e)

    return pos,done,case

def get_case(pos,case):
    if case in pos.__dict__.keys():
        ret = getattr(pos,case)

    else:
        if case=='waist': #Beamwaist as a function of z (z=coordinate along beamline)
            z = pos.zoff
            zR = np.pi*(pos.aim.data.param.w0_laser**2)/pos.aim.data.param.labda
            ret = pos.aim.data.param.w0_laser*((1+((z/zR)**2))**0.5)

        elif case=='coor_starttele': #Done
            ret=pos.coor_starttele__sun

        elif case=='vec_startbeam__send':
            ret = pos.coor_startbeam__send[0]

        elif case=='vec_startbeam__sun':
            v = pos.coor_startbeam__send[0]
            ret = pos.aim.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t_start,v,reverse=False)

        elif case=='vec_startbeam__rec':
            v = pos.vec_startbeam__sun
            ret = pos.aim.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t_end,v,reverse=True)
        
        elif case=='aberration_effect':
            ret = LA.angle(pos.vec_startbeam__rec,pos.vec_startbeam__send)

        elif case=='coor_end':
            ret = pos.coor_endtele__sun

        elif case=='start':
            ret = pos.aim.data.putp(pos.i_send,pos.t_start)+pos.aim.data.param.L_tele*pos.coor_starttele__sun[0]

        elif case=='end':
            ret = pos.aim.data.putp(pos.i_rec,pos.t_end)+pos.aim.data.param.L_tele*pos.coor_endtele__sun[0]

        elif case=='startend__send':
            ret = pos.arm__send

        elif case=='startend__rec':
            v = pos.startend__sun
            ret = pos.aim.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t_end,v,reverse=True)

        elif case=='r':
            ret = (pos.xoff**2+pos.yoff**2)**0.5

        elif case=='R_vec_tele_send':
            ret = LA.matmul(pos.coor_starttele,pos.R_vec_beam__send)

        elif case=='alpha':
            ret = (pos.angx_wf_rec**2+pos.angy_wf_rec**2)**0.5

        elif case=='Ival':
            ret = (pos.I0*np.exp((-2*(pos.xoff**2+pos.yoff**2))/(pos.waist**2)))*(pos.aim.data.param.w0_laser/pos.waist)

        elif case=='I':
            if pos.alpha>pos.aim.data.param.FOV/2.0:
                ret=0.0
            else:
                ret = pos.Ival

        elif case=='I0':
            ret = (pos.aim.data.param.P_L*np.pi*(pos.aim.data.param.w0_laser**2))/4.0

        elif case=='tdel':
            ret = pos.t_end-pos.t_start

        elif case=='invside':
            if pos.side=='l':
                ret = 'r'
            elif pos.side=='r':
                ret='l'

        elif case=='side_send':
            if pos.mode=='send':
                ret = pos.side
            elif pos.mode=='rec':
                ret = pos.invside

        elif case=='side_rec':
            if pos.mode=='send':
                ret = pos.invside
            elif pos.mode=='rec':
                ret = pos.side

    try:
        if pos.t_start>=pos.aim.data.t_all[2] and pos.t_end<=pos.aim.data.t_all[-2]: #...to coarse
            setattr(pos,case,ret)
        else:
            setattr(pos,case,ret*np.nan)
    except UnboundLocalError, e:
        raise UnboundLocalError(e)

    return pos
