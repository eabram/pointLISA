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
        
        elif case=='precision': #...adjust
            ret = 0 

        elif case=='R': #The radius of curvasture R as a function of z
            if pos.precision==0:
                ret = lambda z: pos.zoff
            elif pos.precision>0:
                #z = pos.piston
                if pos.precision==1:
                    ret = lambda z: z
                elif pos.precision==2:
                    #if z!=np.nan:
                    zR = np.pi*(pos.aim.data.param.w0_laser**2)/pos.aim.data.param.labda
                    ret  = lambda z: abs(z*(1+((zR/z)**2)))
                    #else:
                    #    ret = np.nan

        elif case=='end_ksi':
            ret = pos.ksi[0]*pos.coor_end[2]+pos.ksi[1]*pos.coor_end[1]

        elif case=='offset': #Done
            if pos.side=='l':
                ret = pos.offset_l
            elif pos.side=='r':
                ret = pos.offset_r

        elif case=='off': #Done
            ret = LA.matmul(pos.coor_startbeam__send,pos.arm__send)
        elif case=='xoff': #Done
            ret = pos.off[2]
        elif case=='yoff': #Done
            ret = pos.off[1]
        elif case=='zoff': #Done
            ret = pos.off[0]
     
        elif case=='tele_ang': #Done
            if pos.side=='l':
                ret = pos.aim.tele_l_ang(pos.i_self,pos.t)
            elif pos.side=='r':
                ret = pos.aim.tele_r_ang(pos.i_self,pos.t)

        elif case=='PAAM_ang': #Done
            if pos.side=='l':
                ret = pos.aim.beam_l_ang(pos.i_self,pos.t)
            elif pos.side=='r':
                ret = pos.aim.beam_r_ang(pos.i_self,pos.t)
        
        elif case=='invside': #Done
            if pos.side=='l':
                ret='r'
            elif pos.side=='r':
                ret='l'

        elif case=='i_opp': #Done
            if pos.side=='l':
                ret=pos.i_left
            elif pos.side=='r':
                ret=pos.i_right

        elif case=='coor_starttele':
            ret=calc.get_coor_tele(pos.aim,pos.i_send,pos.t_start,pos.side_send,tele_angle=pos.tele_angle_start)

        elif case=='coor_startbeam__send':
            ret=calc.get_coor_beam_out__send(pos.aim,pos.i_send,pos.t_start,pos.side_send,tele_angle=pos.tele_angle_start,beam_angle=pos.beam_angle_start,offset=pos.offset_start)

        elif case=='vec_startbeam__send':
            ret = pos.coor_startbeam__send[0]

        elif case=='vec_startbeam__sun':
            v = pos.coor_startbeam__send
            ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t_end,v,reverse=False)

        elif case=='vec_startbeam__rec':
            v = pos.coor_startbeam__sun
            ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t_end,v,reverse=True)
        
        elif case=='coor_endbeam__send':
            ret = calc.get_coor_beam_out__send(pos.aim,pos.i_send,pos.t_end,pos.side_send,tele_angle=pos.tele_angle_start,beam_angle=pos.beam_angle_start,offset=pos.offset_start)

        elif case=='vec_endbeam__send': #send frame, Sun CS
            ret = pos.coor_endbeam__send[0]

        elif case=='vec_endbeam__sun': #Sun frame, Sun CS
            v = pos.coor_endbeam__send
            ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t_end,v,reverse=False)                

        elif case=='vec_endbeam__rec': #Rec frame, Sun CS
            v = pos.coor_startbeam__sun
            ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t_end,v,reverse=True)

        elif case=='coor_end':
            ret = calc.get_coor_tele(pos.aim,pos.i_rec,pos.t_end,pos.side_rec,tele_angle=pos.tele_angle_end)

        elif case=='start':
            ret = calc.get_start_calc(pos.aim,pos.i_send,pos.t_start,pos.side_send,pos.tele_angle_start)

        elif case=='end':
            ret = calc.get_start_calc(pos.aim,pos.i_rec,pos.t_end,pos.side_rec,pos.tele_angle_end)
            if pos.ksi!=[0,0]:
                ret = ret+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]

        elif case=='startend__sun':
            ret = pos.end-pos.start

        elif case=='startend__send':
            v = pos.startend__sun
            ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t_start,v,reverse=True)

        elif case=='startend__rec':
            v = pos.startend__sun
            ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t_end,v,reverse=True)

#...hier gebleven




        elif case=='r':
            ret = (pos.xoff**2+pos.yoff**2)**0.5

        elif case=='arm__rec':
            ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t_end,pos.startend__sun,reverse=True)
       
        elif case=='arm__send':
            ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t_start,pos.startend__sun,reverse=True) #.......check on time it emits!!!

        elif case=='arm_tele_rec':
            ret = LA.matmul(pos.coor_end,-pos.arm__rec)
        
        elif case=='arm_tele_send':
            ret = LA.matmul(pos.coor_starttele,pos.arm__send)

        elif case=='angx_arm_tele_rec':
            ret = np.arctan(pos.arm_tele_rec[2]/pos.arm_tele_rec[0])

        elif case=='angy_arm_tele_rec':
            ret = np.arctan(pos.arm_tele_rec[1]/pos.arm_tele_rec[0])
        
        elif case=='ang_arm_tele_rec':
            ret = np.arctan(((pos.arm_tele_rec[2]**2+pos.arm_tele_rec[1]**2)**0.5)/pos.arm_tele_rec[0])

        elif case=='angx_arm_tele_send':
                    ret = np.arctan(pos.arm_tele_send[2]/pos.arm_tele_send[0])

        elif case=='angy_arm_tele_send':
            ret = np.arctan(pos.arm_tele_send[1]/pos.arm_tele_send[0])





        elif case=='R_vec_beam__send': # in Sun coordinate system and send inertial frame!!!
            vec = np.array([(pos.R(pos.zoff)**2-pos.xoff**2-pos.yoff**2)**0.5,pos.yoff,pos.xoff]) # In beam frame ...adjust: have to use piston
            ret = LA.matmul(np.linalg.inv(pos.coor_startbeam__send),vec) 

        elif case=='R_vec_beam__sun':
            v = pos.R_vec_beam__send
            ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t_start,v,reverse=False)

        elif case=='R_vec_beam__rec':
            v = pos.R_vec_beam__sun
            ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t_end,v,reverse=True)

        elif case=='R_vec_tele_send':
            ret = LA.matmul(pos.coor_starttele,R_vec_beam__send)

        elif case=='R_vec_tele_rec':
            ret = LA.matmul(pos.coor_end,pos.R_vec_beam__rec)

        elif case=='angx_R_vec_tele_rec':
            ret = np.arctan(pos.R_vec_tele_rec[2]/pos.R_vec_tele_rec[0])

        elif case=='angy_R_vec_tele_rec':
            ret = np.arctan(pos.R_vec_tele_rec[1]/pos.R_vec_tele_rec[0])

        elif case=='ang_R_vec_tele_rec':
            ret = np.arctan(((pos.R_vec_tele_rec[2]**2+pos.R_vec_tele_rec[1]**2)**0.5)/pos.R_vec_tele_rec[0])








        elif case=='beam_receive_rec':
            ret = LA.matmul(pos.coor_end,-pos.vec_endbeam__rec)

        elif case=='angx_rec':
            ret = np.sign(pos.beam_receive_rec[2])*abs(np.arctan(pos.beam_receive_rec[2]/pos.beam_receive_rec[0]))

        elif case=='angy_rec':
            ret = np.sign(pos.beam_receive_rec[1])*abs(np.arctan(pos.beam_receive_rec[1]/pos.beam_receive_rec[0]))

        elif case=='beam_receive_send':
            ret = LA.matmul(pos.coor_starttele,pos.vec_endbeam__send)
        
        elif case=='angx_send':
            ret = np.sign(pos.beam_receive_send[2])*abs(np.arctan(pos.beam_receive_send[2]/pos.beam_receive_send[0]))

        elif case=='angy_send':
            ret = np.sign(pos.beam_receive_send[1])*abs(np.arctan(pos.beam_receive_send[1]/pos.beam_receive_send[0]))

      
        elif case=='tele_vec': #This vector is reversed (pointed away from receiving telescope)
            ret = LA.matmul(pos.coor_starttele,-pos.coor_end[0])

        elif case=='R_vec_beam_send__send': # in send coordinate system and send inertial frame!!!
            v_sun = pos.R_vec_beam__send
            v_send = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t_start,v_sun,reverse=True)
            ret = LA.matmul(pos.coor_startbeam__send,v_send)        

        elif case=='angx_R':
            ret = np.sign(pos.R_vec_beam_send__send[2])*abs(np.arctan(pos.R_vec_beam_send__send[2]/pos.R_vec_beam_send__send[0]))

        elif case=='angy_R':
            ret = np.sign(pos.R_vec_beam_send__send[1])*abs(np.arctan(pos.R_vec_beam_send__send[1]/pos.R_vec_beam_send__send[0]))
       
        elif case=='angx_wf_rec':
            ret = np.sign(pos.R_vec_tele_rec[2])*abs(np.arctan(pos.R_vec_tele_rec[2]/pos.R_vec_tele_rec[0]))

        elif case=='angy_wf_rec':
            ret = np.sign(pos.R_vec_tele_rec[1])*abs(np.arctan(pos.R_vec_tele_rec[1]/pos.R_vec_tele_rec[0]))

        elif case=='angx_wf_send':
            endtele = -pos.coor_endtele[0]
            endtele_send = LA.matmul(pos.coor_starttele,endtele)
            R = pos.R_vec_beam_send__send
            angx_tele = np.sign(endtele_send[2])*np.arctan(abs(endtele_send[2]/endtele_send[0]))
            angx_R = np.sign(R[2])*np.arctan(abs(R[2]/R[0]))
            
            ret = angx_R - angx_tele

        elif case=='angy_wf_send':
            endtele = -pos.coor_endtele[0]
            endtele_send = LA.matmul(pos.coor_starttele,endtele)
            R = pos.R_vec_beam_send__send
            angx_tele = np.sign(endtele_send[1])*np.arctan(abs(endtele_send[1]/endtele_send[0]))
            angx_R = np.sign(R[1])*np.arctan(abs(R[1]/R[0]))

            ret = angx_R - angx_tele
         

        elif case=='dz': #...to do:check
            '''Solves the photon traveling distance from the point of transmitting to the point [x,y,z]'''
            try:
                x = pos.xoff
                y = pos.yoff
                z = pos.zoff
                if z!=np.nan:
                    f_solve = lambda dz: (pos.R(z+dz) - (pos.R(z+dz)**2 - (x**2+y**2))**0.5) - dz
                    ret = scipy.optimize.brentq(f_solve,-0.5*z,0.5*z,xtol=1e-64)
                else:
                    ret=np.nan
                    raise ValueError

            except RuntimeError:
                ret = np.nan

        elif case=='piston': #...adjust
            ret = pos.zoff+pos.dz

        elif case=='alpha':
            ret = (abs(pos.angx_wf_rec)**2 +abs(pos.angy_wf_rec)**2)**0.5
        
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

#    def mean_var(self,i,t,side,ret,mode='mean',Nbins=False,tele_angle_l=False,tele_angle_r=False,beam_angle_l=False,beam_angle_r=False):
#        '''Returns the output value for center (value at the center of the aperture), mean (mean value over the aperture), var (variance of the value over the aperture, mean_var (the mean and variance over the aperture) or mean_surface (matrix of values which represents the aperture pixels)'''
#        if Nbins!=False:
#            self.pupil(Nbins=Nbins)
#        else:
#            try:
#                self.xlist
#            except AttributeError:
#                self.pupil()
#        if type(ret)!=list:
#            ret=[ret]
#
#        func = lambda x,y: calc.values(self,i,t,side,ret=ret,ksi=[x,y],tele_angle_l=tele_angle_l,tele_angle_r=tele_angle_r,beam_angle_l=beam_angle_l,beam_angle_r=beam_angle_r)
#
#        if mode=='center':
#            return getattr(func(0,0),ret[0])
#        elif 'var' in mode or 'mean' in mode:
#            A=[]
#            for x in self.xlist:
#                for y in self.ylist:
#                    if ((x**2+y**2)**0.5)>=self.aim.data.param.D/2.0:
#                        A.append(getattr(func(x,y),ret[0])*np.nan)
#                    else:
#                        A.append(getattr(func(x,y),ret[0]))
#            if mode=='mean':
#                return np.nanmean(A)
#            elif mode=='var':
#                return np.nanvar(A)/(np.float64(len(A) - A.count(np.nan)))
#            elif mode=='mean_var':
#                return np.array([np.nanmean(A),np.nanvar(A)/(np.float64(len(A) - A.count(np.nan)))])
#            elif mode=='mean_surface':
#                return A
#
#### Write and save functions/values
#    
#    def t_calc(self,calc=False,**kwargs):
#        '''Time array'''
#        if calc==True:
#            t0= kwargs.pop('t0',False)
#            tend= kwargs.pop('tend',False)
#            dt= kwargs.pop('dt',False)
#            n= kwargs.pop('n',False)
#
#            if dt==False:
#                try:
#                    dt = self.dt
#                except:
#                    dt = self.aim.data.t_all[1]-self.aim.data.t_all[0]
#            if n!=False:
#                tend = dt*n
#
#            if t0==False:
#                t0 = self.aim.data.t_all[3]
#            if tend==False:
#                tend = self.aim.data.t_all[-3]
#            n = int(np.round((tend-t0)/dt))+1
#            t_plot = np.linspace(t0,tend,n)
#        
#        elif calc==False:
#            t_plot = self.aim.data.t_all
#
#        return t_plot
#    
#    def clear_functions(self):
#        try:
#            del self.func
#        except AttributeError:
#            pass
#        try:
#            del self.sampled
#        except AttributeError:
#            pass
#        
#        return 0
#
#
#    def make_functions(self,include=[],exclude=[],option='both',i='all',side=['l','r'],auto_clear=False,t=False,mode='mean_var',**kwargs):
#        '''Obtains the returned properties'''
#
#        Nbins=kwargs.pop('Nbins',False)
#
#        if auto_clear==True:
#            self.clear_functions()
#
#        ret=[]
#        if include=='all' and exclude!='all':
#            for k in OUTPUT.__dict__.keys():
#                if 'get_' in k:
#                    add = k.split('get_')[-1]
#                    if add not in exclude:
#                        ret.append(add)
#        
#        else:
#            for inc in include:
#                ret.append(inc)
#         
#        func=pointLISA.utils.Object() 
#        sampled=pointLISA.utils.Object()
#        
#        if type(t)==bool:
#            if t==False:
#                t_plot = self.t_calc(calc=True) #add parameters
#        else:
#            t_plot=np.array(t)
#        
#
#        if i =='all':
#            i=range(1,4)
#        elif type(i)!=list:
#            i=[i]
#        if type(side)!=list:
#            side=[side]
#
#        for s in side:
#            setattr(sampled,s,utils.Object())
#            for i_sel in i:
#                setattr(getattr(sampled,s),'i'+str(i_sel),utils.Object())
#
#        ret_new = []
#        for k in ret:
#            ret_new.append(k.replace("'",""))
#        ret = ret_new
#        del ret_new
#
#        for k in ret:
#            print(k)
#            
#            if 'adjust' == k:
#                for s in side:
#                    for i_sel in i:
#                        if self.aim.PAAM_deg==1:
#                            if s=='l':
#                                A = [np.array(getattr(self.aim,'t_adjust')[0][int(i_sel)-1])]
#                                A.append(np.array(getattr(self.aim,'tele_adjust')[0][int(i_sel)-1])) 
#                            elif s=='r':
#                                A = [np.array(getattr(self.aim,'t_adjust')[1][int(i_sel)-1])]
#                                A.append(np.array(getattr(self.aim,'tele_adjust')[1][int(i_sel)-1])) 
#                        elif self.aim.PAAM_deg==2:
#                            if s=='l':
#                                A = [np.array(getattr(self.aim,'t_adjust')[0][int(i_sel)-1])]
#                                A.append(np.array(getattr(self.aim,'tele_adjust')[0][int(i_sel)-1])) 
#                            elif s=='r':
#                                A = [np.array(getattr(self.aim,'t_adjust')[1][int(i_sel)-1])]
#                                A.append(np.array(getattr(self.aim,'tele_adjust')[1][int(i_sel)-1]))   
#                        B = [A,'value='+k+', mode='+str(mode)]
#                        setattr(getattr(getattr(sampled,s),'i'+str(i_sel)),k,B)
#            else:
#                if option=='both' or option=='function':
#                    setattr(func,k,lambda i,t,side: self.mean_var(i,t,side,[k],Nbins=Nbins,mode=mode))
#                if option=='both' or option=='sampled':
#                    for s in side:
#                        for i_sel in i:
#                            A=[t_plot]
#                            try:
#                                A.append(np.array([self.mean_var(i_sel,t,s,ret=[k],Nbins=Nbins,mode=mode) for t in t_plot]))
#                            except TypeError,e:
#                                if "'dict' object is not callable" in str(e):
#                                    A = getattr(self.aim,k)[s][i_sel]
#                            B = [A,'value='+k+', mode='+str(mode)]
#                            setattr(getattr(getattr(sampled,s),'i'+str(i_sel)),k,B)
#        return [func,sampled]
