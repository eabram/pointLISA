from pointLISA import *

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
    #print('Getting: '+case)
    if case not in pos.__dict__.keys():
        if case=='waist': #Beamwaist as a function of z (z=coordinate along beamline)
            z = pos.zoff
            zR = np.pi*(pos.aim.data.param.w0_laser**2)/pos.aim.data.param.labda
            ret = pos.aim.data.param.w0_laser*((1+((z/zR)**2))**0.5)
        
        #elif case=='end_ksi':
        #    ret = pos.ksi[0]*pos.coor_end[2]+pos.ksi[1]*pos.coor_end[1]

#        elif case=='offset_send': #Done
#            ret = pos.offset_send
#        elif case=='offset_rec': #Done
#            ret = pos.offset_rec
#        elif case=='off': #Done
#            ret = pos.off
#        elif case=='xoff': #Done
#            ret = pos.off[2]
#        elif case=='yoff': #Done
#            ret = pos.off[1]
#        elif case=='zoff': #Done
#            ret = pos.off[0]
#     
#        elif case=='tele_ang': #Done
#            if pos.side=='l':
#                ret = pos.aim.tele_l_ang(pos.i_self,pos.t)
#            elif pos.side=='r':
#                ret = pos.aim.tele_r_ang(pos.i_self,pos.t)
#
#        elif case=='PAAM_ang': #Done
#            if pos.side=='l':
#                ret = pos.aim.beam_l_ang(pos.i_self,pos.t)
#            elif pos.side=='r':
#                ret = pos.aim.beam_r_ang(pos.i_self,pos.t)
#        
#        elif case=='invside': #Done
#            if pos.side=='l':
#                ret='r'
#            elif pos.side=='r':
#                ret='l'
#
#        elif case=='i_opp': #Done
#            if pos.side=='l':
#                ret=pos.i_left
#            elif pos.side=='r':
#                ret=pos.i_right

        elif case=='coor_starttele': #Done
            ret=pos.coor_starttele__sun

        elif case=='vec_startbeam__send':
            ret = pos.coor_startbeam__send[0]

        elif case=='vec_startbeam__sun':
            v = pos.coor_startbeam__send[0]
            ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t_start,v,reverse=False)

        elif case=='vec_startbeam__rec':
            v = pos.vec_startbeam__sun
            ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t_end,v,reverse=True)
        
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
            ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t_end,v,reverse=True)

#...hier gebleven




        elif case=='r':
            ret = (pos.xoff**2+pos.yoff**2)**0.5

        elif case=='R_vec_tele_send':
            ret = LA.matmul(pos.coor_starttele,R_vec_beam__send)

        elif case=='alpha':
            (pos.angx_wf_rec**2+pos.angy_wf_rec**2)**0.5

        elif case=='FOVlim':
            if pos.alpha>pos.aim.data.param.FOV:#pos.waist>pos.aim.data.FOV:
                ret=0
            else:
                ret=1

        elif case=='Ival':
            ret = (pos.I0*np.exp((-2*(pos.xoff**2+pos.yoff**2))/(pos.waist**2)))*(pos.aim.data.param.w0_laser/pos.waist)

        elif case=='Ivalx':
            ret = (pos.I0*np.exp((-2*(pos.xoff**2))/(pos.waist**2)))*(pos.aim.data.param.w0_laser/pos.waist)
     
        elif case=='I':
            ret = pos.Ival*pos.FOVlim

        elif case=='I0':
            ret = (pos.aim.data.param.P_L*np.pi*(pos.aim.data.param.w0_laser**2))/2.0

        elif case=='power': #...adjust
            try:
                xlist=pos.xlist
                ylist=pos.ylist
            except AttributeError:
                pos.pupil(Nbins=pos.Nbins)
                xlist=pos.xlist
                ylist=pos.ylist

            if len(xlist)==1 and len(ylist)==1:
                dksi = (self.D**2)*(np.pi/4.0)
            else:
                dksi = (xlist[1]-xlist[0])*(ylist[1]-ylist[0])
            ret = dksi*pos.I

        elif case=='tdel':
            if pos.mode=='send':
                if pos.side=='l':
                    ret = pos.aim.data.L_sl_func_tot(pos.i_self,pos.t)
                elif pos.side=='r':
                    ret = pos.aim.data.L_sr_func_tot(pos.i_self,pos.t)
            elif pos.mode=='rec':
                if pos.side=='l':
                    ret = pos.aim.data.L_rl_func_tot(pos.i_self,pos,t)
                elif pos.side=='r':
                    ret = pos.aim.data.L_rr_func_tot(pos.i_self,pos,t)

        elif case=='tdel0':
            if pos.calc_method=='Abram':
                ret=0
            elif pos.calc_method=='Waluschka':
                ret = pos.tdel

#    elif case=='start':
#        if pos.mode=='send':
#            ret = pos.t
#        elif pos.mode=='rec':
#            ret = pos.t-pos.tdel
#    
#    elif case=='t_end':
#        if pos.mode=='send':
#            ret = pos.t+pos.tdel
#        elif pos.mode=='rec':
#            ret = pos.t

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

        elif case=='tele_angle_start':
            if pos.side_send=='l':
                if pos.getter is True:
                    ret = pos.aim.tele_l_ang(pos.i_send,pos.t_start)
                else:
                    ret = pos.tele_angle_l
            if pos.side_send=='r':
                if pos.getter is True:
                    ret = pos.aim.tele_r_ang(pos.i_send,pos.t_start)
                else:
                    ret = pos.tele_angle_r

        elif case=='tele_angle_end':
            if pos.side_rec=='l':
                if pos.getter is True:
                    ret = pos.aim.tele_l_ang(pos.i_rec,pos.t_end)
                else:
                    ret = pos.tele_angle_l
            if pos.side_rec=='r':
                if pos.getter is True:
                    ret = pos.aim.tele_r_ang(pos.i_rec,pos.t_end)
                else:
                    ret = pos.tele_angle_r

        elif case=='beam_angle_start':
            if pos.side_send=='l':
                if pos.getter is True:
                    ret = pos.aim.beam_l_ang(pos.i_send,pos.t_start)
                else:
                    ret = pos.beam_angle_l
            if pos.side_send=='r':
                if pos.getter is True:
                    ret = pos.aim.beam_r_ang(pos.i_send,pos.t_start)
                else:
                    ret = pos.beam_angle_r

        elif case=='beam_angle_end':
            if pos.side_rec=='l':
                if pos.getter is True:
                    ret = pos.aim.beam_l_ang(pos.i_rec,pos.t_end)
                else:
                    ret = pos.beam_angle_l 
            if pos.side_rec=='r':
                if pos.getter is True:
                    ret = pos.aim.beam_r_ang(pos.i_rec,pos.t_end)
                else:
                    ret = pos.beam_angle_r

        elif case=='offset_start':
            ret = calc.get_offset(pos.aim,pos.i_send,pos.t_start,pos.side_send)

        elif case=='offset_end':
            ret = calc.get_offset(pos.aim,pos.i_end,pos.t_end,pos.side_rec)

        elif case=='coor_endtele':
            if pos.mode=='send':
                if pos.side=='l':
                    ret = calc.get_coor_tele(pos.aim,pos.i_left,pos.t+pos.tdel,'r',tele_angle=pos.tele_angle_r)
                elif pos.side=='r':
                    ret = calc.get_coor_tele(pos.aim,pos.i_right,pos.t+pos.tdel,'l',tele_angle=pos.tele_angle_l)
            elif pos.mode=='rec':
                if pos.side=='l':
                    ret = calc.get_coor_tele(pos.aim,pos.i_self,pos.t-pos.tdel0,'l',tele_angle=pos.tele_angle_l)
                elif pos.side=='r':
                    ret = calc.get_coor_tele(pos.aim,pos.i_self,pos.t-pos.tdel0,'r',tele_angle=pos.tele_angle_r)

        elif case=='t_start':
            if pos.mode=='send':
                ret = pos.t
            elif pos.mode=='rec':
                ret = pos.t-pos.tdel

        elif case=='t_end':
            if pos.mode=='send':
                ret = pos.t+pos.tdel
            elif pos.mode=='rec':
                ret = pos.t
        try: 
            setattr(pos,case,ret)
        except UnboundLocalError, e:
            print(case)
            raise UnboundLocalError(e)
    return pos
