from pointLISA import * 
# This class obtaines the pointing of the PAAM and telescope

class AIM():
    def __init__(self,data=None,setting=utils.Object(),filename=False,**kwargs):        
        # Get settings and parameters
        for k in CALC.__dict__.keys():
            if '__'!=k[0:2]:
                setattr(self,k,getattr(CALC(),k))
        if data!=None:
            for k in kwargs.keys():
                if k in setting.__dict__.keys():
                    delattr(setting,k)
                    setattr(setting,k,kwargs[k])
            self.aimset = setting
            if self.aimset.aim_object!=None:
                self.aimset.tele_control='AIM_object'
                self.aimset.PAAM_control='AIM_object'
                self.aimset.offset_tele='AIM_object'
            self.data = data
            
            if data!=False:
                if self.aimset.import_file==None:
                    print('Start calculating telescope and PAAM aim')
                    self.get_offset_inplane()
            else:
                pass #...adjust
        
                

    def get_offset_inplane(self,**kwargs):
        '''This function obtains the offset between the inplane telescope alignment and the transmitting beam (offset)'''
        print('Getting initial offset angles')
        if self.aimset.offset_tele=='AIM_object':
            self.offset = self.aimset.aim_object.offset
        elif self.aimset.PAAM_deg==1:
            if self.aimset.offset_tele==False:
                offset = lambda i,t,side: 0
            elif self.aimset.offset_tele==True:
                offset = lambda i,t,side: (2*np.pi*self.data.param.L_arm)/(c*year2sec)
            self.offset = offset

        return 0

################# TELESCOPE POINTING #################
    def get_selections(self,i,t,side,mode='send'):
        [i_self,i_left,i_right] = utils.const.i_slr(i)
        if mode=='send':
            t_start = t
            i_send = i_self
            if side=='l':
                i_rec = i_left
                t_end = t_start + self.data.L_sl(i_send,t_start)
            elif side=='r':
                i_rec = i_right
                t_end = t_start + self.data.L_sr(i_send,t_start)
        elif mode=='rec':
            t_end = t
            i_rec = i_self
            if side=='l':
                i_send = i_left
                t_start = t_end - self.data.L_rl(i_rec,t_end)
            elif side=='r':
                i_send = i_right
                t_start = t_end - self.data.L_rr(i_rec,t_end)
        
        return [i_send,i_rec,t_start,t_end,mode]

    def get_output(self,i,t,side,mode='send',tele_l=np.radians(-30.0),tele_r=np.radians(30.0),beam_l=0.0,beam_r=0.0,solve=False,ret='xoff',tele_SS_l=False,tele_SS_r=False):
        [i_send,i_rec,t_start,t_end,mode] = self.get_selections(i,t,side,mode)
        
        if solve is False:
            if (side=='l' and mode=='send') or (side=='r' and mode=='rec'):
                if tele_SS_l!=False:
                    tele_l = tele_SS_l
                else:
                    tele_l = self.tele_l_ang(i_send,t_start)
                if tele_SS_r!=False:
                    tele_r = tele_SS_r
                else:
                    tele_r = self.tele_r_ang(i_rec,t_end)
                beam_l = self.beam_l_ang(i_send,t_start)
                beam_r = self.beam_r_ang(i_rec,t_end)
            elif (side=='r' and mode=='send') or (side=='l' and mode=='rec'):
                tele_l = self.tele_l_ang(i_rec,t_end)
                tele_r = self.tele_r_ang(i_send,t_start)
                beam_l = self.beam_l_ang(i_rec,t_end)
                beam_r = self.beam_r_ang(i_send,t_start)
               
        if (side=='l' and mode=='send') or (side=='r' and mode=='rec'):
            tele_send = tele_l
            tele_rec = tele_r
            beam_send = beam_l
            beam_rec = beam_r
            offset_send = self.offset(i_send,t_start,'l')
            offset_rec = self.offset(i_rec,t_end,'r')
        elif (side=='r' and mode=='send') or (side=='l' and mode=='rec'):
            tele_send = tele_r
            tele_rec = tele_l
            beam_send = beam_r
            beam_rec = beam_l
            offset_send = self.offset(i_send,t_start,'r')
            offset_rec = self.offset(i_rec,t_end,'l')
        
        try:
            return locals()[ret]
        except KeyError:
            coor_startbeam__send = self.beam_coor_out__send(self.data,i_send,t_start,tele_send,beam_send,offset_send)
            coor_starttele__sun = self.coor_tele(self.data,i_send,t_start,tele_send)
            coor_endtele__sun = self.coor_tele(self.data,i_rec,t_end,tele_rec)
            start = coor_starttele__sun[0]*self.data.param.L_tele+self.data.putp(i_send,t_start)
            end = coor_endtele__sun[0]*self.data.param.L_tele+self.data.putp(i_rec,t_end)
            startend__sun = end - start
            arm__send = self.aberration_beam_coor(self.data,i_send,t_start,startend__sun,reverse=True)
            off = LA.matmul(coor_startbeam__send,arm__send)
            if ret=='off':
                return off
            elif ret=='xoff':
                return off[2]
            elif ret=='yoff':
                return off[1]
            elif ret=='zoff':
                return off[0]
            elif ret=='r':
                return (off[1]**2+off[2]**2)**0.5
            elif 'Ival' in ret or ret=='waist':
                zR = np.pi*(self.data.param.w0_laser**2)/self.data.param.labda
                waist = self.data.param.w0_laser*((1+((off[0]/zR)**2))**0.5)
                if ret=='waist':
                    return waist
                else:
                    I0 = (self.data.param.P_L*np.pi*(self.data.param.w0_laser**2))/2.0
                    if ret=='Ival':
                        return (I0*np.exp((-2*(off[2]**2+off[1]**2))/(waist**2)))*(self.data.param.w0_laser/waist)
                    elif ret=='Ivalx':
                        return (I0*np.exp((-2*(off[2]**2))/(waist**2)))*(self.data.param.w0_laser/waist)
        
            else:
                zoff = off[0]
                yoff = off[1]
                xoff = off[2]
                
                zR = (np.pi*(self.data.param.w0_laser**2))/self.data.param.labda
                R_func = lambda dz: (zoff+dz)*(1+((zR/(zoff+dz))**2))
                f_solve = lambda dz: R_func(dz) - ((R_func(dz)**2-xoff**2-yoff**2)**0.5)- dz
                dz = scipy.optimize.brentq(f_solve,-zoff*0.1,zoff*0.1,xtol=1e-64)
                piston = zoff+dz
                R = R_func(dz)
                vec = np.array([(R**2-xoff**2-yoff**2)**0.5,yoff,xoff])
                coor_startbeam__send = self.beam_coor_out__send(self.data,i_send,t_start,tele_send,beam_send,offset_send)
                R_vec_beam__send = LA.matmul(np.linalg.inv(coor_startbeam__send),vec)
                R_vec_beam__sun = self.aberration_beam_coor(self.data,i_send,t_start,R_vec_beam__send,reverse=False)
                R_vec_beam__rec = self.aberration_beam_coor(self.data,i_rec,t_end,R_vec_beam__sun,reverse=True)
                R_vec_tele_rec = LA.matmul(coor_endtele__sun,R_vec_beam__rec)
                if ret=='R_vec_tele_rec':
                    return R_vec_tele_rec
                if ret=='angx_wf_rec':
                    return np.sign(R_vec_tele_rec[2])*abs(np.arctan(R_vec_tele_rec[2]/R_vec_tele_rec[0]))
                elif ret=='angy_wf_rec':
                    return np.sign(R_vec_tele_rec[1])*abs(np.arctan(R_vec_tele_rec[1]/R_vec_tele_rec[0]))
                elif ret=='alpha':
                    return abs(np.arctan(((R_vec_tele_rec[1]**2+R_vec_tele_rec[2]**2)**0.5)/R_vec_tele_rec[0]))
                elif ret=='aberration_effect':
                    return LA.angle(R_vec_beam__send,R_vec_beam__rec)
                elif ret=='piston':
                    return piston
                elif ret=='R':
                    return R
                else:
                    raise ValueError('Please select a proper output parameter')
    
    def get_beam_angle(self,i,t,side,beam_l0 = 0.0,beam_r0 = 0.0,conv_lim=1e-9,loop=1,option=None):
        if option==None:
            option = self.aimset.option_PAAM

        lim = np.radians(20.0)
        [i_send,i_rec,t_start,t_end,mode] = self.get_selections(i,t,side,'send')
        beam_l_new = [beam_l0]
        beam_r_new = [beam_r0]
        conv = [1.0]
        done = False
        l = 0
        while done is False or l<loop:
            if side=='l':
                try:
                    tele_l
                except NameError:
                    tele_l = self.tele_l_ang(i_send,t_start)
                    tele_r = self.tele_r_ang(i_rec,t_end)
                if option=='center':
                    pos_send = lambda beam_l: self.get_output(i_send,t_start,'l',tele_l = tele_l,tele_r=tele_r,beam_l=beam_l,beam_r = 0.0,solve=True,ret='yoff')
                elif option=='wavefront':
                    pos_send = lambda beam_l: self.get_output(i_send,t_start,'l',tele_l = tele_l,tele_r=tele_r,beam_l=beam_l,beam_r = 0.0,solve=True,ret='angy_wf_rec')
                beam_l_new.append(scipy.optimize.brentq(pos_send,-lim+beam_l0,lim+beam_l0))
                conv.append(abs(beam_l_new[-1]-beam_l_new[-2]))
            elif side=='r':
                try:
                    tele_l
                except NameError:
                    tele_l = self.tele_l_ang(i_rec,t_end)
                    tele_r = self.tele_r_ang(i_send,t_start)
                if option=='center':
                    pos_send = lambda beam_r: self.get_output(i_send,t_start,'r',tele_l=tele_l,tele_r=tele_r,beam_l=0.0,beam_r = beam_r,solve=True,ret='yoff')
                elif option=='wavefront':
                    pos_send = lambda beam_r: self.get_output(i_send,t_start,'r',tele_l=tele_l,tele_r=tele_r,beam_l=0.0,beam_r = beam_r,solve=True,ret='angy_wf_rec')

                beam_r_new.append(scipy.optimize.brentq(pos_send,-lim+beam_r0,lim+beam_r0))
                conv.append(abs(beam_r_new[-1]-beam_r_new[-2]))
            if conv <=conv_lim or (conv[-1]-conv[-2])/conv[-2]<0.01:
                done = True
            l = l+1

        if side=='l':
            ret = beam_l_new[-1]
        elif side=='r':
            ret = beam_r_new[-1]

        return ret

    def get_tele_SS(self,solve_for='both',i_all=[1,2,3]):
        def get_SS_func(x,y,x_check):
            '''Returns the SS function'''
            A = [t for t in x if t<x_check]
            val = y[len(A)-1]
            return np.float64(val)

        if solve_for=='power': #...Ivalx voor SS werkt alleen bij center (als wavefrot dan Ivaly ook lage waarde (dus yoff hoog)
            ret = ['Ivalx']
            lim = [self.data.param.I_min]
        elif solve_for=='wavefrontangle':
            ret=['angx_wf_rec']
            lim = [self.aimset.FOV]
        elif solve_for=='both':
            ret=['Ivalx','angx_wf_rec']
            lim = [self.data.param.I_min,self.aimset.FOV]

        t0 = self.data.t_all[0]
        tstop = self.data.t_all[2]-10.0 #...adjust: self.data.t_all[-1] -10.0

        tele_l_SS = {}
        tele_r_SS = {}
        for i in i_all:
            t = [t0]
            tele_l=[]
            tele_r=[]
            run=True
            while run:
                print(t[-1]/tstop)
                [i_send,i_rec,t_start,t_end,mode] = self.get_selections(i,t[-1],'l','send')

                tele_l.append(self.get_tele_angle(i_send,t_start,'l',loop=1))
                tele_r.append(self.get_tele_angle(i_rec,t_end,'r',loop=1))

                if t[-1]>=tstop:
                    run=False
                    break
                pos_send = lambda t: abs(self.get_output(i_send,t,'l','send',tele_l = tele_l[-1],tele_r=tele_r[-1],beam_l=False,beam_r=False,solve=True,ret=ret[0])) - lim[0]
                pos_rec = lambda t: abs(self.get_output(i_rec,t,'r','send',tele_l = tele_l[-1],tele_r=tele_r[-1],beam_l=False,beam_r=False,solve=True,ret=ret[0])) - lim[0]
                if len(ret)==2:
                    pos_send2 = lambda t: abs(self.get_output(i_send,t,'l','send',tele_l = tele_l[-1],tele_r=tele_r[-1],beam_l=False,beam_r=False,solve=True,ret=ret[1])) - lim[1]
                    pos_rec2 = lambda t: abs(self.get_output(i_rec,t,'r','send',tele_l = tele_l[-1],tele_r=tele_r[-1],beam_l=False,beam_r=False,solve=True,ret=ret[1])) - lim[1]

                dt=3600
                solved=False
                while solved is False:
                
                    try:
                        A = scipy.optimize.brentq(pos_send,t[-1]+1.0,t[-1]+dt)
                    except ValueError,e:
                        A = np.nan
                    if len(ret)==2:
                        try:
                            C = scipy.optimize.brentq(pos_send2,t[-1]+1.0,t[-1]+dt)
                        except ValueError,e:
                            C = np.nan
                    else:
                        C = np.nan
                    try:
                        B = scipy.optimize.brentq(pos_rec,t[-1]+1.0,t[-1]+dt)
                    except ValueError,e:
                        B = np.nan
                    if len(ret)==2:
                        try:
                            D = scipy.optimize.brentq(pos_rec2,t[-1]+1.0,t[-1]+dt)
                        except ValueError,e:
                            D = np.nan
                    else:
                        D = np.nan
                    
                    if np.isnan([A,B,C,D]).sum()!=4:
                        t_new = np.nanmin([A,B,C,D])
                        solved=True
                    else:
                        dt=dt+3600
                t.append(t_new)
            
            if t[-1]!=tstop:
                t.append(tstop)
                tele_l.append(np.nan)
                tele_r.append(np.nan)

            t_reset = t

            tele_l[0] = np.nan
            tele_r[0] = np.nan
            tele_l_SS[i_send] = [t_reset,tele_l]
            tele_r_SS[i_rec] = [t_reset,tele_r]
        ang_l = lambda i,t: get_SS_func(tele_l_SS[i][0],tele_l_SS[i][1],t)
        ang_r = lambda i,t: get_SS_func(tele_r_SS[i][0],tele_r_SS[i][1],t)

        return [ang_l,ang_r,tele_l_SS,tele_r_SS]

    def get_tele_angle(self,i,t,side,tele_l0 = np.radians(-30.0),tele_r0 = np.radians(30.0),conv_lim=1e-9,loop=1,option=None):
        if option==None:
            option = self.aimset.option_tele

        lim = np.radians(20.0)
        [i_send,i_rec,t_start,t_end,mode] = self.get_selections(i,t,side,'send')
        tele_l_new = [tele_l0]
        tele_r_new = [tele_r0]
        conv = [1.0]
        done = False
        l = 0
        while done is False or l<loop:
            if side=='l':
                if option=='center':
                    pos_send = lambda tele_l: self.get_output(i_send,t_start,'l',tele_l=tele_l,tele_r = tele_r_new[-1],solve=True)
                    tele_l_new.append(scipy.optimize.brentq(pos_send,-lim+tele_l0,lim+tele_l0))
                    pos_rec = lambda tele_r: self.get_output(i_rec,t_end,'r',tele_l=tele_l_new[-1],tele_r = tele_r,solve=True)
                    tele_r_new.append(scipy.optimize.brentq(pos_rec,-lim+tele_r0,lim+tele_r0))
                elif option=='wavefront':
                    pos_send = lambda tele_r: self.get_output(i_send,t_start,'l',tele_l=tele_l_new[-1],tele_r = tele_r,solve=True,ret='angx_wf_rec')
                    tele_r_new.append(scipy.optimize.brentq(pos_send,-lim+tele_r0,lim+tele_r0))
                    pos_rec = lambda tele_l: self.get_output(i_rec,t_end,'r',tele_l=tele_l,tele_r = tele_r_new[-1],solve=True,ret='angx_wf_rec')
                    tele_l_new.append(scipy.optimize.brentq(pos_rec,-lim+tele_l0,lim+tele_l0))
                conv.append(max(abs(tele_l_new[-1]-tele_l_new[-2]),abs(tele_r_new[-1]-tele_r_new[-2])))
            elif side=='r':
                if option=='center':
                    pos_send = lambda tele_r: self.get_output(i_send,t_start,'r',tele_l=tele_l_new[-1],tele_r = tele_r,solve=True)
                    tele_r_new.append(scipy.optimize.brentq(pos_send,-lim+tele_r0,lim+tele_r0))
                    pos_rec = lambda tele_l: self.get_output(i_rec,t_end,'l',tele_l=tele_l,tele_r = tele_r_new[-1],solve=True)
                    tele_l_new.append(scipy.optimize.brentq(pos_rec,-lim+tele_l0,lim+tele_l0))
                elif option=='wavefront':
                    pos_send = lambda tele_l: self.get_output(i_send,t_start,'r',tele_l=tele_l,tele_r = tele_r_new[-1],solve=True,ret='angx_wf_rec')
                    tele_l_new.append(scipy.optimize.brentq(pos_send,-lim+tele_l0,lim+tele_l0))
                    pos_rec = lambda tele_r: self.get_output(i_rec,t_end,'l',tele_l=tele_l_new[-1],tele_r = tele_r,solve=True,ret='angx_wf_rec')
                    tele_r_new.append(scipy.optimize.brentq(pos_rec,-lim+tele_r0,lim+tele_r0))
                conv.append(max(abs(tele_l_new[-1]-tele_l_new[-2]),abs(tele_r_new[-1]-tele_r_new[-2])))
            #print(conv[-1])
            if conv <=conv_lim or (conv[-1]-conv[-2])/conv[-2]<0.01:
                done = True
            l = l+1

        if side=='l':
            #print(pos_rec(tele_l_new[-1]))
            ret = tele_l_new[-1]
        elif side=='r':
            #print(pos_rec(tele_r_new[-1]))
            ret = tele_r_new[-1]

        return ret

    def tele_control_ang_fc(self,option=None,value=False,lim=False):
        '''Obtains the telescope pointing angles for a continuous actuation (full_control)'''
        # Option 'wavefront' means poiting with the purpose of getting a zero/small tilt of the receiving wavefront
        # 'center' means pointing it to te center of the receiving telescope aperture
        if option==None:
            option = self.aimset.option_tele

        print('Telescope pointing strategy: '+option)
        
        method = self.aimset.solve_method
        if method=='fast' and option=='center':
            ang_l = lambda i,t: - LA.angle(self.data.v_l_in(i,t),self.data.r_func(i,t)) - self.offset(i,t,'l')
            ang_r = lambda i,t: LA.angle(self.data.v_r_in(i,t),self.data.r_func(i,t)) - self.offset(i,t,'r')

        else:
            max_count=1 #...adjust for better optimization
            
            ang_l = lambda i,t: self.get_tele_angle(i,t,'l',loop=max_count)
            ang_r = lambda i,t: self.get_tele_angle(i,t,'r',loop=max_count)
 
        self.aimset.option_tele = option
        self.aimset.tele_control = 'full_control'

        return [ang_l,ang_r]


    def tele_aim(self,lim=1e-10):
        '''Obtains the telescope pointing angles (for the selected telescope pointing method)'''
        method=self.aimset.tele_control
        option = self.aimset.option_tele

        try:
            print('The telescope control method is: '+method)
        except:
            print('The telescope control method is: user defined')

        print(' ')

        if self.aimset.import_file==None:
            if method=='no_control':
                # For no_control (no pointing)
                self.tele_l_ang_func = lambda i,t: np.radians(-30)
                self.tele_r_ang_func = lambda i,t: np.radians(30)

            elif method=='full_control':
                # For ful_control (continuous poiting)
                [self.tele_ang_l_fc,self.tele_ang_r_fc] = self.tele_control_ang_fc()
                self.tele_l_ang_func = self.tele_ang_l_fc
                self.tele_r_ang_func = self.tele_ang_r_fc

            elif 'SS' in method:
                # For Step-and_Stare
                [self.tele_l_ang_func,self.tele_r_ang_func,self.adjust_l,self.adjust_r] = self.get_tele_SS()

            elif type(method)==list and method[0]=='Imported pointing': #...to do
                # Import (previously obtained) angles
                print(method[0])
                self.tele_l_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'l',x=method[1]['SC'+str(i)+', left']['x'],y=method[1]['SC'+str(i)+', left']['y'])
                self.tele_r_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'r',x=method[1]['SC'+str(i)+', right']['x'],y=method[1]['SC'+str(i)+', right']['y'])

                t_adjust={}
                tele_ang_adjust = {}
                for i in range(1,4):
                    t_adjust[str(i)]={}
                    tele_ang_adjust[str(i)]={}
                    t_adjust[str(i)]['l']=method[1]['SC'+str(i)+', left']['x']
                    t_adjust[str(i)]['r']=method[1]['SC'+str(i)+', right']['x']
                    tele_ang_adjust[str(i)]['l']=method[1]['SC'+str(i)+', left']['y']
                    tele_ang_adjust[str(i)]['r']=method[1]['SC'+str(i)+', right']['y']


                self.t_adjust = t_adjust
                self.tele_ang_adjust = tele_ang_adjust

            elif type(self.aimset.tele_control)==tuple: # and self.tele_control[0].options.tele_control=='SS':
                self.tele_l_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'l',x=getattr(self.aimset.tele_control[0].l,'i'+str(i)).adjust[0][0],y=getattr(self.aimset.tele_control[0].l,'i'+str(i)).adjust[0][1])
                self.tele_r_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'r',x=getattr(self.aimset.tele_control[0].r,'i'+str(i)).adjust[0][0],y=getattr(self.aimset.tele_control[0].r,'i'+str(i)).adjust[0][1])
                
                print('Read imported telescope SS angles')
                t_adjust=[[0,0,0],[0,0,0]]
                tele_ang_adjust = [[0,0,0],[0,0,0]]
                fl = [0,0,0]
                fr = [0,0,0]
                for i in range(1,4):
                    t_adjust[0][i-1]=getattr(self.aimset.tele_control[0].l,'i'+str(i)).adjust[0][0]
                    t_adjust[1][i-1]=getattr(self.aimset.tele_control[0].r,'i'+str(i)).adjust[0][0]
                    tele_ang_adjust[0][i-1]=getattr(self.aimset.tele_control[0].l,'i'+str(i)).adjust[0][1]
                    tele_ang_adjust[1][i-1]=getattr(self.aimset.tele_control[0].r,'i'+str(i)).adjust[0][1]

                self.t_adjust = t_adjust
                self.tele_ang_adjust = tele_ang_adjust
                self.tele_adjust = tele_ang_adjust

            elif self.aimset.tele_control=='AIM_object':
                self.tele_l_ang_func = self.aimset.aim_object.tele_l_ang
                self.tele_r_ang_func = self.aimset.aim_object.tele_r_ang

            else:
                raise ValueError('Please select valid telescope pointing method')
    
        else:
            print('Importing pointing angles from:')
            print(self.aimset.import_file)
            ret = pointLISA.read_write.read_output(filenames=self.aimset.import_file)
            tele_l_ang=[]
            tele_r_ang=[]
            if method=='SS':
                self.tele_l_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'l',x=ret['SC'+str(i)+', left']['x'],y=ret['SC'+str(i)+', left']['y'])
                self.tele_r_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'r',x=ret['SC'+str(i)+', right']['x'],y=ret['SC'+str(i)+', right']['y'])

                t_adjust={}
                tele_ang_adjust = {}
                for i in range(1,4):
                    t_adjust[str(i)]={}
                    tele_ang_adjust[str(i)]={}
                    
                    [t_adjust[str(i)]['l'],tele_ang_adjust[str(i)]['l']]=getattr(ret[0].l,'i'+str(i)).adjust[0]
                    [t_adjust[str(i)]['r'],tele_ang_adjust[str(i)]['r']]=getattr(ret[0].r,'i'+str(i)).adjust[0]

                self.tele_l_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'l',x=t_adjust[str(i)]['l'],y=tele_ang_adjust[str(i)]['l'])
                self.tele_r_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'l',x=t_adjust[str(i)]['r'],y=tele_ang_adjust[str(i)]['r'])
                
                self.t_adjust = t_adjust
                self.tele_ang_adjust = tele_ang_adjust

            else:

                for i in range(1,4):
                    t_l =  getattr(getattr(ret[0],'l'),'i'+str(i)).tele_ang
                    t_r =  getattr(getattr(ret[0],'r'),'i'+str(i)).tele_ang
                    tele_l_ang.append(calc.interpolate(t_l[0][0],t_l[0][1]))
                    tele_r_ang.append(calc.interpolate(t_r[0][0],t_r[0][1]))

                self.tele_l_ang = lambda i,t: tele_l_ang[i-1](t)
                self.tele_r_ang = lambda i,t: tele_r_ang[i-1](t)
            
            if self.aimset.PAAM_deg==2:
                try:
                    delattr(self,'offset')
                except AttributeError:
                    pass
                offset={}
                offset['l']={}
                offset['r']={}
                for i in range(1,4):
                    offset_l =  getattr(getattr(ret[0],'l'),'i'+str(i)).offset
                    offset_r =  getattr(getattr(ret[0],'r'),'i'+str(i)).offset
                    offset['l'][i] = calc.interpolate(offset_l[0][0],offset_l[0][1])
                    offset['r'][i] = calc.interpolate(offset_r[0][0],offset_r[0][1])
                self.offset={}
                self.offset = offset

            else:
                self.get_offset_inplane(self.aimset.offset_tele)

        try:
            self.tele_l_ang
        except AttributeError:
            if self.aimset.sampled==True:
                try:
                    self.t_sample
                except AttributeError:
                    self.t_sample = calc.get_t_sample(self,speed=self.aimset.sample_speed)
                
                tele_l_ang=[]
                tele_r_ang=[]
                print("Sampling and fitting telescope angles") #...add offset voor PAAM_deg==2
                for i in range(1,4):
                    tele_l_ang.append(calc.interpolate(self.t_sample['l'][i-1],np.array([self.tele_l_ang_func(i,t) for t in self.t_sample['l'][i-1]])))
                    tele_r_ang.append(calc.interpolate(self.t_sample['r'][i-1],np.array([self.tele_r_ang_func(i,t) for t in self.t_sample['r'][i-1]])))
                self.tele_l_ang_samp = lambda i,t: tele_l_ang[i-1](t)
                self.tele_r_ang_samp = lambda i,t: tele_r_ang[i-1](t)
                
                self.tele_l_ang = self.tele_l_ang_samp
                self.tele_r_ang = self.tele_r_ang_samp

            else:
                self.tele_l_ang = self.tele_l_ang_func
                self.tele_r_ang = self.tele_r_ang_func

        return 0























    def PAAM_control_ang_fc(self,option=None,tele_l_ang=False,tele_r_ang=False):
        '''Obtains the PAAM pointing angles for a continuous actuation (full_control)'''
        if option==None:
            option = self.aimset.option_PAAM
        print('PAAM pointing strategy: '+option)

        method = self.aimset.solve_method
        if method=='fast' and option=='center':
            ang_l = lambda i,t: -0.5*self.data.PAA.out(i,t,'l')
            ang_r = lambda i,t: 0.5*self.data.PAA.out(i,t,'r')

        max_count=1 #...adjust for better optimization

        ang_l = lambda i,t: self.get_beam_angle(i,t,'l',loop=max_count)
        ang_r = lambda i,t: self.get_beam_angle(i,t,'r',loop=max_count)

        self.aimset.option_PAAM = option
        self.aimset.PAAM_control = 'full_control'

        return [ang_l,ang_r]


    def PAAM_aim(self):
        '''Obtains the PAAM pointing angles (for the selected telescope pointing method)'''
        method=self.aimset.PAAM_control
        option = self.aimset.option_PAAM
        print('The PAAM control method is: ' +method)
        print(' ')

        # Obtaining PAAM angles for 'fc' (full_control), 'nc' (no_control) and 'SS' (step and stare)
        
        if self.aimset.import_file==None:
            if method=='full_control':
                [ang_l,ang_r] = self.PAAM_control_ang_fc(option=option)

            elif method=='no_control':
                ang_l = lambda i,t: 0
                ang_r = lambda i,t: 0

            elif method=='SS':
                ang_l_SS = lambda i,t: ang_fc_l(i,t-(t%dt)) # Adjusting the pointing every dt seconds
                ang_r_SS = lambda i,t: ang_fc_r(i,t-(t%dt))
                print('Taken '+method+' step response for PAAM SS control with tau='+str(tau)+' sec')
                mode='overdamped'
            elif method=='AIM_object':
                ang_l=self.aimset.aim_object.beam_l_ang
                ang_r=self.aimset.aim_object.beam_r_ang

            self.beam_l_ang = ang_l
            self.beam_r_ang = ang_r
                
        else:
            print('Importing pointing angles from:')
            print(self.aimset.import_file)
            ret = pointLISA.read_write.read_output(filenames=self.data.import_file)
            beam_l_ang=[]
            beam_r_ang=[]
            for i in range(1,4):
                b_l =  getattr(getattr(ret[0],'l'),'i'+str(i)).PAAM_ang
                b_r =  getattr(getattr(ret[0],'r'),'i'+str(i)).PAAM_ang
                beam_l_ang.append(calc.interpolate(b_l[0][0],b_l[0][1]))
                beam_r_ang.append(calc.interpolate(b_r[0][0],b_r[0][1]))

            self.beam_l_ang_func = lambda i,t: beam_l_ang[i-1](t)
            self.beam_r_ang_func = lambda i,t: beam_r_ang[i-1](t)
            self.beam_l_ang = lambda i,t: beam_l_ang[i-1](t)
            self.beam_r_ang = lambda i,t: beam_r_ang[i-1](t)

        return self
 

    def twoPAAM_fc(self,aim_dummy):
        aim_dummy.tele_aim()
        aim_dummy.PAAM_aim()

        add_l = lambda i,t: aim_dummy.get_output(i,t,'l',mode='rec',ret='angx_wf_rec')
        add_r = lambda i,t: aim_dummy.get_output(i,t,'r',mode='rec',ret='angx_wf_rec')

        offset={}
        tele_l_ang = lambda i,t: aim_dummy.tele_l_ang(i,t)-add_l(i,t)
        offset['l'] = lambda i,t: add_l(i,t) +aim_dummy.offset(i,t,'l')
        tele_r_ang = lambda i,t: aim_dummy.tele_r_ang(i,t)-add_r(i,t)
        offset['r'] = lambda i,t: add_r(i,t)+aim_dummy.offset(i,t,'r')
        beam_l_ang = aim_dummy.beam_l_ang
        beam_r_ang = aim_dummy.beam_r_ang
        #offset=lambda i,t,s: offset[s](i,t)
        return [tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang,offset]

    def twoPAAM(self):
        ### On;y works with option_tele==center and option_PAAM==center
        kwargs = {'PAAM_deg':1,'option_tele':'center','option_PAAM':'center','tele_control':'full_control','PAAM_control':'full_control'}
        kwargs['solve_method'] = self.aimset.solve_method
        aimset_dummy = const.get_settings(settings_input=self.data.settings,select='aimset',kwargs=kwargs)
        aim_dummy = AIM(input_file=None,data=self.data,setting = aimset_dummy,filename=False)
        self.aim_dummy = aim_dummy

        [tele_l_ang_fc,tele_r_ang_fc,beam_l_ang_fc,beam_r_ang_fc,offset_fc] = self.twoPAAM_fc(aim_dummy)
        offset={}
        if self.aimset.tele_control=='no_control':
            self.tele_l_ang = lambda i,t: np.radians(-30.0) 
            self.tele_r_ang = lambda i,t: np.radians(30.0)
            offset['l'] = lambda i,t: aim_dummy.tele_l_ang(i,t) - self.tele_l_ang(i,t) + aim_dummy.offset(i,t,'l')
            offset['r'] = lambda i,t: aim_dummy.tele_r_ang(i,t) - self.tele_r_ang(i,t) + aim_dummy.offset(i,t,'r')
            self.offset = lambda i,t,s: offset[s](i,t)

        elif self.aimset.tele_control=='full_control':
            self.tele_l_ang = tele_l_ang_fc
            self.tele_r_ang = tele_r_ang_fc
            self.offset = lambda i,t,s: offset_fc[s](i,t)

        elif self.aimset.tele_control=='SS':
            pass #...adjust

        if self.aimset.PAAM_control=='no_control':
            self.beam_l_ang = lambda i,t: 0.0
            self.beam_r_ang = lambda i,t: 0.0

        elif self.aimset.PAAM_control=='full_control':
            self.beam_l_ang = beam_l_ang_fc
            self.beam_r_ang = beam_r_ang_fc

        return 0

class CALC():
    # Changes of coordinate system
    def coor_SC(self,data,i,t):
        '''Returns the coordinates of a spacecraft in [r,n,x] components'''
        t_calc=t

        r = LA.unit(data.r_func(i,t_calc))
        n = LA.unit(data.n_func(i,t_calc))
        x = np.cross(n,r)

        return np.array([r,n,x])

    def coor_tele(self,data,i,t,ang_tele):
        '''Returns the coordinate system of telescope (same as SC but rotated over ang_tele inplane)'''
        L_tele = data.param.L_tele
        [r,n,x] = self.coor_SC(data,i,t)
        tele = r*L_tele
        tele = LA.rotate(tele,n,ang_tele)
        r = LA.unit(tele)
        x = np.cross(n,r)

        return np.array([r,n,x])

    def aberration_beam_coor(self,data,i,t,v,reverse=False): # if reverse==True: SUN-->SC, if reverse==False: SC-->SUN
        if data.stat.aberration==False:
            ret = v
        elif data.stat.aberration==True:
            V = data.vel.abs(i,t)
            if reverse==True:
                V=-V
            v_mag = np.linalg.norm(v)
            c_vec = LA.unit(v)*data.param.c
            ret = LA.unit(c_vec+V)*v_mag

        return ret

    def beam_coor_out__send(self,data,i,t,ang_tele,ang_paam,ang_tele_offset): # beam coordinates as seen from send frame, Sun coordinate
        '''Retunrs the coordinate system of the transmitted beam (same as SC but rotated over ang_tele inplane and ang_tele outplane)'''
        [r,n,x] = self.coor_tele(data,i,t,ang_tele+ang_tele_offset) #Telescope coordinate system

        r = LA.unit(LA.rotate(r,x,ang_paam)) # Rotate r in out of plane over ang_paam
        n = np.cross(r,x)

        return np.array([r,n,x])

