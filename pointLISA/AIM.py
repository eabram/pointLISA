from pointLISA import * 
# This class obtaines the pointing of the PAAM and telescope

class AIM():
    def __init__(self,data=None,setting=utils.Object(),filename=False,**kwargs):        
        # Get settings and parameters
        if data!=None:
            self.aimset = setting
            self.data = data
            
            if data!=False:
                if self.aimset.import_file==None:
                    print('Start calculating telescope and PAAM aim')
                    self.get_offset_inplane(self.aimset.offset_tele)
            else:
                pass #...adjust
                

    def get_offset_inplane(self,option,**kwargs):
        '''This function obtains the offset between the inplane telescope alignment and the transmitting beam (offset)'''
        print('Getting initial offset angles')
        if self.aimset.PAAM_deg==1:
            offset = {'l': {1: 0.0, 2: 0.0, 3: 0.0},
     'r': {1: 0.0, 2: 0.0, 3: 0.0}}
            
            if option=='0':
                option =0

            if option==0:
                pass

            elif option==True:
                offset = {'l': {1: 0.0, 2: 0.0, 3: 0.0},
     'r': {1: 0.0, 2: 0.0, 3: 0.0}}

            elif '/' in option:
                ret = calc.read(direct=option)
                for k1 in ret:
                    for k2 in ret[k1]:
                        for k3 in ret[k1][k2]:
                            for k4 in ret[k1][k2][k3]:
                                inp = ret[k1][k2][k3][k4]['angx_func_rec mean']
                                break

                for side in ['l','r']:
                    for SC in range(1,4):
                        if side=='l':
                            label='SC'+str(SC)+', left'
                        elif side=='r':
                            label='SC'+str(SC)+', right'
                        ang = inp[label]['y'][3:-3]

                        offset[side][SC] = -np.nanmean(ang)
            
            elif 'read'==option:
                if 'filename' in kwargs.keys():
                    read_file = kwargs['filename']

                else:
                    read_folder = os.path.dirname(os.path.realpath(__file__))+'/parameters/'+self.data.stat.calc_method+'/'

                    print(read_folder)
                    for (dirpath, dirnames, filenames) in os.walk(read_folder):
                        for f in filenames:
                            read_file = open(dirpath+f)
                            break
                    out=''
                    count=0
                    for line in read_file:
                        count=count+1
                        if count>1:
                            out=out+line
                    out.replace('\n','')
                    offset_all = yaml.load(out)
                    keys = offset_all.keys()
                    keys.sort()
                    days = len(self.data.t_all)
                    try:
                        offset = offset_all[days]
                    except KeyError:
                        try:
                            pos_l = calc.nearest_smaller_value(keys,days)
                            k_l = keys[pos_l]
                            k_r = keys[pos_l+1]
                            sm = offset_all[k_l]
                            bi = offset_all[k_r]
                            offset={}
                            for s in sm.keys():
                                offset[s]={}
                                for i in sm[s].keys():
                                    if days==k_l:
                                        offset[s][i] = sm[s][i]
                                    else:
                                        offset[s][i] = sm[s][i] + ((bi[s][i]-sm[s][i])/(k_r-k_l))*(days-k_l)
                        except:
                            offset = offset_all[keys[-1]]

                    print('Offset is:')
                    print(offset)
                    print('')
                    read_file.close()
            
            elif type(option==dict):
                offset = option

            else:
                raise ValueError("Please select offset tele values or method")

        elif self.aimset.PAAM_deg==2:
            print('Two axis PAAM selected')
            offset={'l': {1: 0.0, 2: 0.0, 3: 0.0},
     'r': {1: 0.0, 2: 0.0, 3: 0.0}}
            self.offset_init=True

        self.offset = offset

        return 0



################# TELESCOPE POINTING #################

#    def tele_point_calc(self,i,t,side,option,lim=False,method=False,value=0,scale=1,max_count=20,tele_l0=None,tele_r0=None,beam_l0=None,beam_r0=None,offset_l0=None,offset_r0=None,**kwargs):
#        '''Calculates the (full control) telescope pointing angles (with the center or wavefront method)'''
#        [i_self,i_left,i_right] = const.i_slr(i)
#        if option=='center':
#            if lim==False:
#                lim = self.aimset.limit_xoff
#            if side=='l':
#                ang = self.tele_center_calc(i,t,lim=lim,value=value,tele_l=tele_l0,tele_r=tele_r0,beam_l=beam_l0,beam_r=beam_r0,offset_l=offset_l0,offset_r=offset_r0)[0][0]
#            elif side=='r':
#                ang = self.tele_center_calc(const.i_slr(i)[2],t,lim=lim,value=value,tele_l=tele_l0,tele_r=tele_r0,beam_l=beam_l0,beam_r=beam_r0,offset_l=offset_l0,offset_r=offset_r0)[0][1]
#
#        elif option=='wavefront':
#            if method==False:
#                method = self.aimset.tele_method_solve
#            if lim==False:
#                lim=self.aimset.limit_angx
#
#            if side=='l':
#                ang = self.get_tele_wavefront(i,t,'l',method,scale=scale,lim=lim,max_count=max_count,value=value,tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_l=beam_l0,beam_r=beam_r0,offset_l=offset_l0,offset_r=offset_r0)
#            elif side=='r':
#                ang = self.get_tele_wavefront(i_right,t,'r',method,scale=scale,lim=lim,max_count=max_count,value=value,tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_l=beam_l0,beam_r=beam_r0,offset_l=offset_l0,offset_r=offset_r0)
#
#        return ang
#
#    def tele_center_calc(self,i,t,scale=1,lim=1e-12,max_count=5,value=0,tele_l=False,tele_r=False,beam_l=False,beam_r=False,offset_l=False,offset_r=False):
#        '''Obtains the telescope pointing angle when the telesope is pointed with the center method'''
#        [i_self,i_left,i_right] = const.i_slr(i)
#
#        lim = np.radians(5.0)
#        if tele_l is False:
#            tele_l=self.tele_l_ang(i_self,t)
#        elif tele_l==None:
#            tele_l=np.radians(np.float64(-30.0))
#        if tele_r is False:
#            tele_r=self.tele_r_ang(i_left,t)
#        elif tele_r==None:
#            tele_r=np.radians(np.float64(30.0))
#        if beam_l is False:
#            beam_l=self.beam_l_ang(i_self,t)
#        elif beam_l==None:
#            beam_l=np.float64(0.0)
#        if beam_r is False:
#            beam_r=self.beam_r_ang(i_self,t)
#        elif beam_r==None:
#            beam_r=np.float64(0.0)
#
#        tele_l_old = tele_l
#        tele_r_old = tele_r
#
#        pos_send = lambda tele_l: calc.values(self,i_self,t,'l',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['xoff']).xoff
#        send_solve = lambda tele_l: pos_send(tele_l)-value
#
#        try:
#            tele_l_new = scipy.optimize.brentq(send_solve,-lim-np.radians(30.0),lim-np.radians(30.0))
#        except ValueError,e:
#            if str(e)=='f(a) and f(b) must have different signs':
#                tele_l_new=np.nan
#
#        if tele_l_new!=np.nan:
#            pos_rec = lambda tele_r: calc.values(self,i_left,t,'r',tele_angle_l=tele_l_new,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['xoff']).xoff
#            rec_solve = lambda tele_r: pos_rec(tele_r)-value
#
#            try:
#                tele_r_new = scipy.optimize.brentq(rec_solve,-lim+np.radians(30.0),lim+np.radians(30.0))
#            except ValueError,e:
#                if str(e)=='f(a) and f(b) must have different signs':
#                    tele_r_new=np.nan
#        else:
#            tele_r_new = np.nan
#        
#        #print(pos_send(tele_l_new),pos_rec(tele_r_new))
#        return [[tele_l_new,tele_r_new], False]

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
                t_start = t_end + self.data.L_rr(i_rec,t_end)
        
        return [i_send,i_rec,t_start,t_end,mode]

    def get_output(self,i,t,side,mode='send',tele_l=np.radians(-30.0),tele_r=np.radians(30.0),beam_l=0.0,beam_r=0.0,offset_l=0.0,offset_r=0.0,solve=False,ret='xoff'):
        [i_send,i_rec,t_start,t_end,mode] = self.get_selections(i,t,side,mode)
        
        if solve is False:
            if (side=='l' and mode=='send') or (side=='r' and mode=='rec'):
                tele_l = self.tele_l_ang(i_send,t_start)
                tele_r = self.tele_r_ang(i_rec,t_end)
                beam_l = self.beam_l_ang(i_send,t_start)
                beam_r = self.beam_r_ang(i_rec,t_end)
                offset_l = self.offset['l'][i_send] #...adjust
                offset_r = self.offset['r'][i_rec] #...adjust
            elif (side=='r' and mode=='send') or (side=='l' and mode=='rec'):
                tele_l = self.tele_l_ang(i_rec,t_end)
                tele_r = self.tele_r_ang(i_send,t_start)
                beam_l = self.beam_l_ang(i_rec,t_end)
                beam_r = self.beam_r_ang(i_send,t_start)
                offset_l = self.offset['l'][i_rec] #...adjust
                offset_r = self.offset['r'][i_send] #...adjust
            
        if side=='l':
            tele_send = tele_l
            tele_rec = tele_r
            beam_send = beam_l
            beam_rec = beam_r
            offset_send = offset_l
            offset_rec = offset_r
        elif side=='r':
            tele_send = tele_r
            tele_rec = tele_l
            beam_send = beam_r
            beam_rec = beam_l
            offset_send = offset_r
            offset_rec = offset_l
            
        coor_startbeam__send = calc.beam_coor_out__send(self.data,i_send,t_start,tele_send,beam_send,offset_send)
        coor_starttele__sun = calc.coor_tele(self.data,i_send,t_start,tele_send)
        coor_endtele__sun = calc.coor_tele(self.data,i_rec,t_end,tele_rec)
        start = coor_starttele__sun[0]*self.data.param.L_tele+self.data.putp(i_send,t_start)
        end = coor_endtele__sun[0]*self.data.param.L_tele+self.data.putp(i_rec,t_end)
        startend__sun = end - start
        arm__send = calc.aberration_beam_coor(self.data,i_send,t_start,startend__sun,reverse=True)
        off = LA.matmul(coor_startbeam__send,arm__send)
        if ret=='off':
            return off
        elif ret=='xoff':
            return off[2]
        elif ret=='yoff':
            return off[1]
        elif ret=='zoff':
            return off[0]
    
        else:
            [zoff,yoff,xoff] = off
            zR = np.pi*(self.data.param.w0_laser**2)/self.data.param.labda
            R = abs(zoff*(1+((zR/zoff)**2)))
            vec = np.array([(R**2-xoff**2-yoff**2)**0.5,yoff,xoff]) # In beam frame ...adjust: have to use piston
            coor_startbeam__send = calc.beam_coor_out__send(self.data,i_send,t_start,tele_send,beam_send,offset_send)
            R_vec_beam__send = LA.matmul(np.linalg.inv(coor_startbeam__send),vec)
            R_vec_beam__sun = calc.aberration_beam_coor(self.data,i_send,t_start,R_vec_beam__send,reverse=False)
            R_vec_beam__rec = calc.aberration_beam_coor(self.data,i_rec,t_end,R_vec_beam__sun,reverse=True)
            R_vec_tele_rec = LA.matmul(coor_endtele__sun,R_vec_beam__rec)
            if ret=='R_vec_tele_rec':
                return R_vec_tele_rec
            if ret=='angx_wf_rec':
                return np.sign(R_vec_tele_rec[2])*abs(np.arctan(R_vec_tele_rec[2]/R_vec_tele_rec[0]))
            elif ret=='angy_wf_rec':
                return np.sign(R_vec_tele_rec[1])*abs(np.arctan(R_vec_tele_rec[1]/R_vec_tele_rec[0]))
            else:
                raise ValueError('Please select a proper output parameter')


    def get_tele_center(self,i,t,side,tele_l0 = np.radians(-30.0),tele_r0 = np.radians(30.0),conv_lim=1e-9,loop=1):
        lim = np.radians(5)   
        [i_send,i_rec,t_start,t_end,mode] = self.get_selections(i,t,side,'send')
        tele_l_new = [tele_l0]
        tele_r_new = [tele_r0]
        conv = [1.0]
        done = False
        l = 0
        while done is False or l<loop:
            pos_send = lambda tele_l: self.get_output(i_send,t_start,'l',tele_l=tele_l,tele_r = tele_r_new[-1],solve=True)
            tele_l_new.append(scipy.optimize.brentq(pos_send,-lim+tele_l0,lim+tele_l0))
            pos_rec = lambda tele_r: self.get_output(i_rec,t_end,'r',tele_l=tele_l_new[-1],tele_r = tele_r,solve=True)
            tele_r_new.append(scipy.optimize.brentq(pos_rec,-lim+tele_r0,lim+tele_r0))
            conv.append(max(abs(tele_l_new[-1]-tele_l_new[-2]),abs(tele_r_new[-1]-tele_r_new[-2])))
            #print(conv[-1])
            if conv <=conv_lim or (conv[-1]-conv[-2])/conv[-2]<0.01:
                done = True
            l = l+1
                
        if side=='l':
            ret = tele_l_new[-1]
        elif side=='r':
            ret = tele_r_new[-1]
        
        return ret

    def get_tele_wavefront(self,i,t,side,tele_l0 = np.radians(-30.0),tele_r0 = np.radians(30.0),conv_lim=1e-9,loop=1):
        lim = np.radians(20.0)
        [i_send,i_rec,t_start,t_end,mode] = self.get_selections(i,t,side,'send')
        tele_l_new = [tele_l0]
        tele_r_new = [tele_r0]
        conv = [1.0]
        done = False
        l = 0
        while done is False or l<loop:
            if side=='l':
                pos_send = lambda tele_r: self.get_output(i_send,t_start,'l',tele_l=tele_l_new[-1],tele_r = tele_r,solve=True,ret='angx_wf_rec')
                tele_r_new.append(scipy.optimize.brentq(pos_send,-lim+tele_r0,lim+tele_r0))
                pos_rec = lambda tele_l: self.get_output(i_rec,t_end,'r',tele_l=tele_l,tele_r = tele_r_new[-1],solve=True,ret='angx_wf_rec')
                tele_l_new.append(scipy.optimize.brentq(pos_rec,-lim+tele_l0,lim+tele_l0))
                conv.append(max(abs(tele_l_new[-1]-tele_l_new[-2]),abs(tele_r_new[-1]-tele_r_new[-2])))
            elif side=='r':
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




#    def get_tele_wavefront(self,i,t,side,method,scale=1,lim=1e-12,max_count=20,print_on=False,value=0.0,tele_angle_l=None,tele_angle_r=None,beam_l=None,beam_r=None,offset_l=False,offset_r=False):
#        '''Gets all telescope pointing angles along one arm for the wavefront method'''
#        if side=='l':
#            i_l = i
#            tdel=0
#        elif side=='r':
#            i_l = const.i_slr(i)[2]
#            i_r = i
#            tdel = self.data.L_rr(i_r,t)
#
#        tele_l_old = 0.0
#        tele_r_old = 0.0
#        if tele_angle_l==None:
#            tele_l = np.radians(-30.0)
#        else:
#            tele_l = tele_angle_l
#        if tele_angle_r==None:
#            tele_r = np.radians(30.0)
#        else:
#            tele_r = tele_angle_r
#
#        count=0
#        while count<max_count:
#            [[tele_l_new,tele_r_new],con] = self.tele_wavefront_calc(i,t,tele_l=tele_l,tele_r=tele_r,beam_l=beam_l,beam_r=beam_r,offset_l=offset_l,offset_r=offset_r)
#            count = count+1
#            if count>= max_count:
#                mode = 'Maximum iteration limit has been reached'
#                tele_l = tele_l_new
#                tele_r = tele_r_new
#                if print_on:
#                    print(mode)
#                break
#            elif max(con)<1.0e-9: #max(tele_l_new-tele_l,tele_r_new-tele_r)<1.0e9:
#                mode = 'Result is converged'
#                tele_l = tele_l_new
#                tele_r = tele_r_new
#                if print_on:
#                    print(mode)
#                break
#            tele_l = tele_l_new
#            tele_r = tele_r_new
#
#        if side=='l':
#            return tele_l
#
#        elif side=='r':
#            return tele_r
#
#    def tele_wavefront_calc(self,i_l,t,scale=1,lim=1e-12,max_count=5,print_on=False,value=0,tele_l=False,tele_r=False,beam_l=False,beam_r=False,offset_l=False,offset_r=False):
#        '''Obtains the telescope pointing angle when the telesope is pointed with the center method'''
#        [i_self,i_left,i_right] = const.i_slr(i_l)
#
#        lim = np.radians(5.0)
#        if tele_l is False:
#            tele_l=self.tele_l_ang(i_self,t)
#        elif tele_l==None:
#            tele_l=np.radians(np.float64(-30.0))
#        if tele_r is False:
#            tele_r=self.tele_r_ang(i_left,t)
#        elif tele_r==None:
#            tele_r=np.radians(np.float64(30.0))
#        if beam_l is False:
#            beam_l=self.beam_l_ang(i_self,t)
#        elif beam_l==None:
#            beam_l=np.float64(0.0)
#        if beam_r is False:
#            beam_r=self.beam_r_ang(i_self,t)
#        elif beam_r==None:
#            beam_r=np.float64(0.0)
#
#        tele_l_old = tele_l
#        tele_r_old = tele_r
#
#        par = 'angx_wf_rec'
#
#        pos_l = getattr(calc.values(self,i_self,t,'l',mode='rec',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=[par]),par)
#        tele_l_new = tele_l-pos_l
#
#        pos_r = getattr(calc.values(self,i_left,t,'r',mode='rec',tele_angle_l=tele_l_new,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=[par]),par)
#        tele_r_new = tele_r-pos_r
#
#        return [[tele_l_new,tele_r_new], [pos_l,pos_r]]




    def tele_control_ang_fc(self,option=None,value=False,lim=False):
        '''Obtains the telescope pointing angles for a continuous actuation (full_control)'''
        # Option 'wavefront' means poiting with the purpose of getting a zero/small tilt of the receiving wavefront
        # 'center' means pointing it to te center of the receiving telescope aperture
        if option==None:
            option = self.aimset.option_tele

        print('Telescope pointing strategy: '+option)
        
        max_count=1
        
        if option=='center':
            ang_l = lambda i,t: self.get_tele_center(i,t,'l',loop=max_count)
            ang_r = lambda i,t: self.get_tele_center(i,t,'r',loop=max_count)
        elif option=='wavefront':
            ang_l = lambda i,t: self.get_tele_wavefront(i,t,'l',loop=max_count)
            ang_r = lambda i,t: self.get_tele_wavefront(i,t,'r',loop=max_count)

        #ang_l = lambda i,t: self.tele_point_calc(i,t,'l',option,max_count=max_count,scale=scale,value=value,tele_l0=np.radians(-30.0),tele_r0=np.radians(30.0),beam_l0=0.0,beam_r0=0.0,lim=lim)
        #ang_r = lambda i,t: self.tele_point_calc(i,t,'r',option,max_count=max_count,scale=scale,value=value,tele_l0=np.radians(-30.0),tele_r0=np.radians(30.0),beam_l0=0.0,beam_r0=0.0,lim=lim)
        
        self.aimset.option_tele = option
        self.aimset.tele_control = 'full_control'

        return [ang_l,ang_r]


    def tele_aim(self,method=False,lim=1e-10,tele_ang_extra=True,option=False):
        '''Obtains the telescope pointing angles (for the selected telescope pointing method)'''
        if method==False:
            method=self.aimset.tele_control
        else:
            self.aimset.tele_control=method
            print('Changed tele_control option')
        
        
        if option==False:
            option = self.aimset.option_tele
        else:
            self.aimset.option_tele = option
            print('Changed option_tele option')

        if option=='wavefront':
            value=self.aimset.value_wavefront
        elif option=='center':
            value=self.aimset.value_center

        try:
            print('The telescope control method is: '+method)
        except:
            print('The telescope control method is: user defined')

        print(' ')

        if self.aimset.import_file==None:
            if method=='no_control':
                # For no_control (no pointing)
                print('tele_ang_extra: '+str(tele_ang_extra))
                if tele_ang_extra==False:
                    offset_l = {'1':0,'2':0,'3':0}
                    offset_r = {'1':0,'2':0,'3':0}
                else:
                    [offset_l,offset_r] = tele_ang_extra

                self.tele_l_ang_func = lambda i,t: np.radians(-30)+offset_l[str(i)]
                self.tele_r_ang_func = lambda i,t: np.radians(30)+offset_r[str(i)]

            elif method=='full_control':
                # For ful_control (continuous poiting)
                [self.tele_ang_l_fc,self.tele_ang_r_fc] = self.tele_control_ang_fc(option=option,value=value)
                self.tele_l_ang_func = self.tele_ang_l_fc
                self.tele_r_ang_func = self.tele_ang_r_fc

            elif 'SS' in method:
                # For Step-and_Stare
                tele_l_tot=[0,0,0]
                tele_r_tot=[0,0,0]
                t_l_adjust=[0,0,0]
                t_r_adjust=[0,0,0]
                tele_l_adjust=[0,0,0]
                tele_r_adjust=[0,0,0]

                if option=='wavefront':
                    ret = 'angx_wf_rec'
                    lim=self.aimset.FOV #Klopt niet want alleen in inplane i.p.v. totaal, kan ook met I ...adjust
                elif option=='center':
                    #ret = 'angx_wf_rec'
                    #lim=self.aimset.FOV #Klopt niet want alleen in inplane i.p.v. totaal, kan ook met I ...adjust
                    ret = 'Ivalx'
                    lim=self.aimset.power/(((self.data.param.D**2)/4.0)*np.pi)
                
                if self.aimset.testSS==False:
                    i_start=1
                    i_end=4
                else:
                    i_start=self.aimset.testSS
                    i_end = i_start+1

                for link in range(i_start,i_end):
                    t_plot = self.data.t_all[2:-3]
                    t_adjust,[tele_l,tele_r],i_left,i_right = calc.SS_value(self,link,t_plot[0],t_plot[-1],'solve',lim,ret=ret,print_on=False,value=value,tele_l=None,tele_r=None,offset_l=False,offset_r=False)
                    f_l = lambda t: calc.get_SS_func(t_adjust,tele_l,t)
                    f_r = lambda t: calc.get_SS_func(t_adjust,tele_r,t)
                    tele_l_tot[i_left-1] = f_l
                    tele_r_tot[i_right-1] = f_r
                    t_l_adjust[i_left-1] = t_adjust
                    t_r_adjust[i_right-1] = t_adjust
                    tele_l_adjust[i_left-1] = tele_l
                    tele_r_adjust[i_right-1] = tele_r

                tele_l_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'l',x=t_l_adjust[i-1],y=tele_l_adjust[i-1])
                tele_r_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'r',x=t_r_adjust[i-1],y=tele_r_adjust[i-1])

                self.t_l_adjust = lambda i,t: t_l_adjust[i-1]
                self.t_r_adjust = lambda i,t: t_r_adjust[i-1]
                self.tele_l_ang_SS = tele_l_ang
                self.tele_r_ang_SS = tele_r_ang

                self.tele_l_ang_func = self.tele_l_ang_SS
                self.tele_r_ang_func = self.tele_r_ang_SS

                self.t_adjust = [t_l_adjust,t_r_adjust]
                self.tele_adjust = [tele_l_adjust,tele_r_adjust]
                self.tele_adjust_samp=[tele_l_adjust,tele_r_adjust]
                delat=['tele_l_ang','tele_r_ang','tele_l_ang_func','tele_r_ang_func']
                for d in delat:
                    try:
                        delattr(self,d)
                    except AttributeError:
                        pass
                self.tele_l_ang_func = self.tele_l_ang_SS
                self.tele_r_ang_func = self.tele_r_ang_SS



            elif type(method)==list and method[0]=='Imported pointing':
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























    def PAAM_control_ang_fc(self,option='center',tele_l_ang=False,tele_r_ang=False):
        '''Obtains the PAAM pointing angles for a continuous actuation (full_control)'''
        print('PAAM pointing strategy: '+option)

        if tele_l_ang==False:
            tele_l_ang = self.tele_l_ang
        if tele_r_ang==False:
            tele_r_ang = self.tele_r_ang

        if option=='center':
            ang_l = lambda i,t: calc.PAAM_center_calc(self,i,t,tele_l=tele_l_ang(i,t),tele_r=tele_r_ang(const.i_slr(i)[1],t),lim=self.aimset.limit_yoff,method=self.aimset.PAAM_method_solve,para=self.aimset.optimize_PAAM,value=self.aimset.optimize_PAAM_value,margin=self.aimset.optimize_PAAM_margin)[0][0]
            ang_r = lambda i,t: calc.PAAM_center_calc(self,const.i_slr(i)[2],t,tele_l=tele_l_ang(const.i_slr(i)[2],t),tele_r=tele_r_ang(i,t),lim=self.aimset.limit_yoff,method=self.aimset.PAAM_method_solve,para=self.aimset.optimize_PAAM,value=self.aimset.optimize_PAAM_value,margin=self.aimset.optimize_PAAM_margin)[0][1]

        elif option=='wavefront':
            ang_l = lambda i,t: calc.PAAM_wavefront_calc(self,i,t,'l',lim=self.aimset.limit_angy,tele_l=self.tele_l_ang(i,t),tele_r=self.tele_r_ang(const.i_slr(i)[1],t))
            ang_r = lambda i,t: calc.PAAM_wavefront_calc(self,i,t,'r',lim=self.aimset.limit_angy,tele_l=self.tele_l_ang(const.i_slr(i)[2],t),tele_r=self.tele_r_ang(i,t))

        self.aimset.PAAM_control = 'full_control'
        self.aimset.option_PAAM = option
        return [ang_l,ang_r]


    def PAAM_aim(self,method=False,dt=3600*24,tau=1,mode='overdamped',PAAM_ang_extra=False,option=False):
        '''Obtains the PAAM pointing angles (for the selected telescope pointing method)'''
        if method==False:
            method = self.aimset.PAAM_control
        else:
            self.aimset.PAAM_control = method
            print('PAAM_method has changed')

        if option==False:
            option = self.aimset.option_PAAM
        else:
            self.aimset.option_PAAM = option
            print('option_PAAM has changed')


        print('The PAAM control method is: ' +method)
        print(' ')

        # Obtaining PAAM angles for 'fc' (full_control), 'nc' (no_control) and 'SS' (step and stair)
        
        if self.aimset.import_file==None:
            try:
                self.beam_l_ang
                print('Reading imported PAAM angles')
                
            except AttributeError:
                if method=='full_control':
                    [ang_l,ang_r] = self.PAAM_control_ang_fc(option=option)

                elif method=='no_control':
                    #self.do_static_tele_angle('PAAM')
                    if PAAM_ang_extra==False:
                        ang_l = lambda i,t: 0
                        ang_r = lambda i,t: 0
                    else:
                        [offset_l,offset_r] = PAAM_ang_extra
                        ang_l = lambda i,t: offset_l[i-1]*0.5
                        ang_r = lambda i,t: offset_r[i-1]*0.5

                elif method=='SS':
                    ang_l_SS = lambda i,t: ang_fc_l(i,t-(t%dt)) # Adjusting the pointing every dt seconds
                    ang_r_SS = lambda i,t: ang_fc_r(i,t-(t%dt))
                    print('Taken '+method+' step response for PAAM SS control with tau='+str(tau)+' sec')
                    mode='overdamped'

                elif method=='SS_lim':
                    ang_l_SS = lambda i: self.SS_control(ang_fc_l,i,False,dt=False,xlim=False,accuracy=3600,FOV=self.FOV_control,step=False)
                    ang_r_SS = lambda i: self.SS_control(ang_fc_r,i,False,dt=False,xlim=False,accuracy=3600,FOV=self.FOV_control,step=False)
                    print('Taken '+method+' step response for PAAM SS control with tau='+str(tau)+' sec and step limit='+str(self.FOV_control*1e6)+' radians')
                    mode='not_damped' #...damped SS not implemented jet for SS_lim
                else:
                    raise ValueError('Please select a valid PAAM pointing method')


                if 'SS' in method:
                    ang_l = self.step_response(ang_l_SS,'PAAM',dt,tau=tau,mode=mode)
                    ang_r = self.step_response(ang_r_SS,'PAAM',dt,tau=tau,mode=mode)
                    f_noise_l = lambda i,t: (ang_l(i,t)-ang_l_SS(i,t))**2
                    f_noise_r = lambda i,t: (ang_r(i,t)-ang_r_SS(i,t))**2
                    self.PAAM_ang_l_SS = ang_l_SS
                    self.PAAM_ang_r_SS = ang_r_SS
                    self.PAAM_step = dt


                self.beam_l_ang = ang_l
                self.beam_r_ang = ang_r
                
                if self.aimset.sampled==True:
                    sampled = self.sample()
                    self.aim_sampled = sampled
                else:
                    self.aim_sampled=False

            if self.aimset.inp!=False:
                self.aim_sampled = self

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

    def sample(self):
        '''Returns a new AIM object with the pointing angles sampled'''
        self_new = AIM(data=None)
        not_copy = ['tele_l_ang','tele_r_ang','beam_l_ang','beam_r_ang']
        for k,value in self.__dict__.items():
            if k not in not_copy:
                setattr(self_new,k,value)
        self_new.sampled=True
        
        try:
            self.t_sample
        except AttributeError:
            self.t_sample = calc.get_t_sample(self,speed=self.aimset.sample_speed)

        tele_l_ang=[]
        tele_r_ang=[]
        beam_l_ang=[]
        beam_r_ang=[]
        try:
            self.tele_l_ang_samp
            tele_samp=True
        except AttributeError:
            tele_samp=False

        for i in range(1,4):
            if tele_samp==False:
                tele_l_ang.append(calc.interpolate(self.t_sample['l'][i-1],np.array([self.tele_l_ang(i,t) for t in self.t_sample['l'][i-1]])))
                tele_r_ang.append(calc.interpolate(self.t_sample['r'][i-1],np.array([self.tele_r_ang(i,t) for t in self.t_sample['l'][i-1]])))
            beam_l_ang.append(calc.interpolate(self.t_sample['l'][i-1],np.array([self.beam_l_ang(i,t) for t in self.t_sample['l'][i-1]])))
            beam_r_ang.append(calc.interpolate(self.t_sample['r'][i-1],np.array([self.beam_r_ang(i,t) for t in self.t_sample['r'][i-1]])))
        if tele_samp==False:
           self_new.tele_l_ang = lambda i,t: tele_l_ang[i-1](t)
           self_new.tele_r_ang = lambda i,t: tele_r_ang[i-1](t)
        else:
           self_new.tele_l_ang = self.tele_l_ang_samp
           self_new.tele_r_ang = self.tele_r_ang_samp

        self_new.beam_l_ang = lambda i,t: beam_l_ang[i-1](t)
        self_new.beam_r_ang = lambda i,t: beam_r_ang[i-1](t)

        return self_new

    def get_tele_coor(self,i,t,tele_l_ang,tele_r_ang):
        ''' Calculating new pointing vectors and telescope coordinate system '''
        tele_l_coor = calc.coor_tele(self.data,i,t,tele_l_ang(i,t))
        tele_r_coor = calc.coor_tele(self.data,i,t,tele_r_ang(i,t))
        tele_l_vec = LA.unit(tele_l_coor[0])*self.data.param.L_tele
        tele_r_vec = LA.unit(tele_r_coor[0])*self.data.param.L_tele
        tele_l_start = tele_l_vec+np.array(self.data.putp(i,t))
        tele_r_start = tele_r_vec+np.array(self.data.putp(i,t))

        return [[tele_l_coor,tele_r_coor],[tele_l_vec,tele_r_vec],[tele_l_start,tele_r_start]]

    def get_beam_coor(self,i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang,offset): #...only send frame
        ''' Calculating new pointing vectors and beam coordinate system '''
        if offset == False:
            offset = self.aimset.offset_tele

        offset_l = calc.get_offset(self,i,t,'l')
        offset_r = calc.get_offset(self,i,t,'r')

        beam_l_coor = calc.beam_coor_out(self.data,i,t,tele_l_ang(i,t),beam_l_ang(i,t),offset_l)
        beam_r_coor = calc.beam_coor_out(self.data,i,t,tele_r_ang(i,t),beam_r_ang(i,t),offset_r)

        # Calculating the Transmitted beam direction and position of the telescope aperture
        beam_l_direction = beam_l_coor[0]
        beam_r_direction = beam_r_coor[0]
        beam_l_start = self.data.param.L_tele*beam_l_direction+np.array(self.data.putp(i,t)) #...kan weg
        beam_r_start = self.data.param.L_tele*beam_r_direction+np.array(self.data.putp(i,t)) #...kan weg

        return [[beam_l_coor,beam_r_coor],[beam_l_direction,beam_r_direction],[beam_l_start,beam_r_start]]

    def get_coordinate_systems(self,iteration_val=False,option='self'):
        '''Make functions of the coordinate vectors and systems'''
        if iteration_val==False:
            tele_l_ang = self.tele_l_ang
            tele_r_ang = self.tele_r_ang
            beam_l_ang = self.beam_l_ang
            beam_r_ang = self.beam_r_ang
            offset = self.aimset.offset_tele
        else:
            [[tele_l_ang,tele_r_ang],[beam_l_ang,beam_r_ang]] = self.get_functions_from_sampling(self.get_sampled_pointing(option=option))

        self.tele_l_coor = lambda i,t: self.get_tele_coor(i,t,tele_l_ang,tele_r_ang)[0][0]
        self.tele_r_coor = lambda i,t: self.get_tele_coor(i,t,tele_l_ang,tele_r_ang)[0][1]
        self.tele_l_vec = lambda i,t: self.get_tele_coor(i,t,tele_l_ang,tele_r_ang)[1][0]
        self.tele_r_vec = lambda i,t: self.get_tele_coor(i,t,tele_l_ang,tele_r_ang)[1][1]
        self.tele_l_start = lambda i,t: self.get_tele_coor(i,t,tele_l_ang,tele_r_ang)[2][0] 
        self.tele_r_start = lambda i,t: self.get_tele_coor(i,t,tele_l_ang,tele_r_ang)[2][1] 


        self.beam_l_coor = lambda i,t: self.get_beam_coor(i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang,offset)[0][0]
        self.beam_r_coor = lambda i,t: self.get_beam_coor(i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang,offset)[0][1]
        self.beam_l_direction = lambda i,t: self.get_beam_coor(i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang,offset)[1][0]
        self.beam_r_direction = lambda i,t: self.get_beam_coor(i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang,offset)[1][1]
        self.beam_l_start = lambda i,t: self.get_beam_coor(i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang,offset)[2][0] #...kan weg
        self.beam_r_start = lambda i,t: self.get_beam_coor(i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang,offset)[2][1] #...kan weg

        self.tele_l_ang_calc = tele_l_ang
        self.tele_r_ang_calc = tele_r_ang
        self.beam_l_ang_calc = beam_l_ang
        self.beam_r_ang_calc = beam_r_ang

        return 0 
   
    def twoPAAM_calc(self,i,t,side):
        '''Calculates the two PAAM pointing angles for a dual axis PAAM'''
        max_count=5
        scale=1
        value = self.aimset.value_center

        [i_self,i_left,i_right] = const.i_slr(i)
        if side=='l':
            if self.aimset.tele_control=='full_control':
                Dt = self.data.L_sl(i_self,t)
                ang_in_l = calc.tele_point_calc(self,i,t,'l','center',max_count=max_count,scale=scale,value=value,tele_l0=np.radians(-30.0),tele_r0=np.radians(30.0),beam_l0=0.0,beam_r0=0.0,offset_l0=0.0,offset_r0=0.0)
                ang_in_r = calc.tele_point_calc(self,const.i_slr(i)[1],t+Dt,'r','center',max_count=max_count,scale=scale,value=value,tele_l0=np.radians(-30.0),tele_r0=np.radians(30.0),beam_l0=0.0,beam_r0=0.0,offset_l0=0.0,offset_r0=0.0)
            elif self.aimset.tele_control=='no_control':
                ang_in_l = np.radians(-30.0)
                ang_in_r = np.radians(30.0)

            if self.aimset.PAAM_control=='full_control':
                ang_out = calc.PAAM_center_calc(self,i,t,tele_l=ang_in_l,tele_r=ang_in_r,lim=self.aimset.limit_yoff,method=self.aimset.PAAM_method_solve,para=self.aimset.optimize_PAAM,value=self.aimset.optimize_PAAM_value,margin=self.aimset.optimize_PAAM_margin,beam_l=0.0,beam_r=0.0,offset_l=0.0,offset_r=0.0)[0][0]
            elif self.aimset.PAAM_control=='no_control':
                ang_out = 0.0

            return [ang_in_l,ang_out]

        elif side=='r':
            if self.aimset.tele_control=='full_control':
                Dt = self.data.L_sr(i_self,t)
                ang_in_l = calc.tele_point_calc(self,const.i_slr(i)[2],t+Dt,'l','center',max_count=max_count,scale=scale,value=value,tele_l0=np.radians(-30.0),tele_r0=np.radians(30.0),beam_l0=0.0,beam_r0=0.0,offset_l0=0.0,offset_r0=0.0)
                ang_in_r = calc.tele_point_calc(self,i,t,'r','center',max_count=max_count,scale=scale,value=value,tele_l0=np.radians(-30.0),tele_r0=np.radians(30.0),beam_l0=0.0,beam_r0=0.0,offset_l0=0.0,offset_r0=0.0)
            elif self.aimset.tele_control=='no_control':
                ang_in_l=np.radians(-30.0)
                ang_in_r=np.radians(30.0)

            if self.aimset.PAAM_control=='full_control':
                ang_out = calc.PAAM_center_calc(self,const.i_slr(i)[2],t,tele_l=ang_in_l,tele_r=ang_in_r,lim=self.aimset.limit_yoff,method=self.aimset.PAAM_method_solve,para=self.aimset.optimize_PAAM,value=self.aimset.optimize_PAAM_value,margin=self.aimset.optimize_PAAM_margin,beam_l=0.0,beam_r=0.0,offset_l=0.0,offset_r=0.0)[0][1]
            elif self.aimset.PAAM_control=='no_control':
                ang_out = 0.0

            return [ang_in_r,ang_out]
    
    def twoPAAM_angles(self,sampled=None,**kwargs):
        '''Runner function for obtaineing the dual axis PAAM poiting angles and accompanying objects'''
        param = ['tele_control','PAAMin_control','PAAMout_control']
        for p in param:
            try:
                print(kwargs[p])
                delattr(self.aimset,p)
                setattr(self.aimset,p,p_new)
                setattr(self,p,p_new)
            except KeyError:
                pass
            try:
                getattr(self,p)
            except AttributeError:
                setattr(self,p,getattr(self.aimset,p))
                pass

        self.twoPAAM_PAAMout_aim()
        self.beam_l_ang = self.beam_l_ang_fc
        self.beam_r_ang = self.beam_r_ang_fc
        
        if self.aimset.tele_control=='SS':
            self.twoPAAM_tele_aim_SS()

        elif self.aimset.tele_control=='full_control':
            self.tele_l_ang = lambda i,t: self.twoPAAM_tele_aim(i,t,'l')[0]
            self.tele_r_ang = lambda i,t: self.twoPAAM_tele_aim(i,t,'r')[0]
            offset={}
            offset['l'] = lambda i,t: self.twoPAAM_tele_aim(i,t,'l')[1]
            offset['r'] = lambda i,t: self.twoPAAM_tele_aim(i,t,'r')[1]
            self.offset=lambda i,t,s: offset[s](i,t)

        elif type(self.aimset.tele_control)==tuple: # and self.aimset.tele_control[0].options.tele_control=='SS':
            self.tele_l_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'l',x=getattr(self.aimset.tele_control[0].l,'i'+str(i)).adjust[0][0],y=getattr(self.aimset.tele_control[0].l,'i'+str(i)).adjust[0][1])
            self.tele_r_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'r',x=getattr(self.aimset.tele_control[0].r,'i'+str(i)).adjust[0][0],y=getattr(self.aimset.tele_control[0].r,'i'+str(i)).adjust[0][1])

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

            offset_calc={}
            offset_calc['l'] = lambda i,t: self.twoPAAM_PAAMin_aim_SS(i,t,'l')
            offset_calc['r'] = lambda i,t: self.twoPAAM_PAAMin_aim_SS(i,t,'r')
            self.offset = lambda i,t,s: offset_calc[s](i,t)

        if (self.aimset.sampled == True and sampled==None) or sampled==True:
            t_sample = self.data.t_all
            tele_func_dict = {}
            beam_func_dict = {}
            t_func_dict={}

            offset = {}
           
            for s in ['l','r']:
                tele_func_dict[s]={}
                beam_func_dict[s]={}
                t_func_dict[s]={}
                offset[s]={}

            for i in range(1,4):
                tele_l_calc = []
                tele_r_calc = []
                beam_l_calc = []
                beam_r_calc = []
                t_l_calc = []
                t_r_calc = []
                offset_l_calc=[]
                offset_r_calc=[]
                for t in t_sample:
                    
                    tele_l_calc.append(self.tele_l_ang(i,t))
                    tele_r_calc.append(self.tele_r_ang(i,t))
                    beam_l_calc.append(self.beam_l_ang(i,t))
                    beam_r_calc.append(self.beam_l_ang(i,t))
                    t_l_calc.append(t)
                    t_r_calc.append(t)
                    
                    offset_l_calc.append(self.offset(i,t,'l'))
                    offset_r_calc.append(self.offset(i,t,'r'))
                     
                    print(t/t_sample[-1])
                
                tele_func_dict['l'][i] = calc.interpolate(np.array(t_l_calc),np.array(tele_l_calc))
                beam_func_dict['l'][i] = calc.interpolate(np.array(t_l_calc),np.array(beam_l_calc))
                tele_func_dict['r'][i] = calc.interpolate(np.array(t_r_calc),np.array(tele_r_calc))
                beam_func_dict['r'][i] = calc.interpolate(np.array(t_r_calc),np.array(beam_r_calc))
#
                offset['l'][i] = calc.interpolate(np.array(t_l_calc),np.array(offset_l_calc))
                offset['r'][i] = calc.interpolate(np.array(t_r_calc),np.array(offset_r_calc))
            
            tele_l_ang = lambda i,t: tele_func_dict['l'][i](t)
            beam_l_ang = lambda i,t: beam_func_dict['l'][i](t)
            tele_r_ang = lambda i,t: tele_func_dict['r'][i](t)
            beam_r_ang = lambda i,t: beam_func_dict['r'][i](t)
            del_list=['tele_l_ang','tele_r_ang','beam_l_ang','beam_r_ang','offset']
            for l in del_list:
                try:
                    delattr(self,l)
                except AttributeError:
                    pass

            self.beam_l_ang = beam_l_ang
            self.beam_r_ang = beam_r_ang
            self.tele_l_ang = tele_l_ang
            self.tele_r_ang = tele_r_ang
            self.offset = lambda i,t,s: offset[s][i](t)

        return 0
    
    def twoPAAM_tele_aim_SS_calc(self,i,t,s,tele_angle,beam_angle_l=False,beam_angle_r=False):
        tele_r0 = np.radians(30)
        tele_l0 = -np.radians(30.0)
        offset_l0=0.0
        offset_r0=0.0
        [i_self,i_left,i_right] = const.i_slr(i)

        if s=='l':
            ret='off'
            s_inv = 'r'
            tdel = self.data.L_rl(i,t)
            if self.data.stat.calc_method=='Waluschka':
                tdel0 = tdel
            elif self.data.stat.calc_method=='Abram':
                tdel0 = 0.0
            
            if beam_angle_l==False:
                beam_angle_l = self.beam_l_ang(i_self,t)
            if beam_angle_r==False:
                beam_angle_r = self.beam_r_ang(i_left,t-tdel)

            A = getattr(calc.values(self,i_left,t-tdel,s_inv,tele_angle_l=tele_angle,tele_angle_r=tele_r0,beam_angle_l=beam_angle_l,beam_angle_r=beam_angle_r,offset_l=offset_l0,offset_r=offset_r0,mode='send',ret=[ret]),ret)

            offset_r = np.sign(A[2])*abs(np.arctan(A[2]/A[0]))

            ret2 = 'angx_wf_rec'
            B = getattr(calc.values(self,i_self,t-tdel0,s,tele_angle_l=tele_angle,tele_angle_r=tele_r0,beam_angle_l=beam_angle_l,beam_angle_r=beam_angle_r,offset_l=offset_l0,offset_r=offset_r0+offset_r,mode='rec',ret=[ret2]),ret2)
        
        elif s=='r':
            ret='off'
            s_inv = 'l'
            tdel = self.data.L_rr(i,t)
            if self.data.stat.calc_method=='Waluschka':
                tdel0 = tdel
            elif self.data.stat.calc_method=='Abram':
                tdel0 = 0.0

            if beam_angle_l==False:
                beam_angle_l = self.beam_l_ang(i_right,t-tdel)
            if beam_angle_r==False:
                beam_angle_r = self.beam_r_ang(i_self,t)
            A = getattr(calc.values(self,i_right,t-tdel,s_inv,tele_angle_l=tele_l0,tele_angle_r=tele_angle,beam_angle_l=beam_angle_l,beam_angle_r=beam_angle_r,offset_l=offset_l0,offset_r=offset_r0,mode='send',ret=[ret]),ret)

            offset_l = np.sign(A[2])*abs(np.arctan(A[2]/A[0]))

            ret2 = 'angx_wf_rec'
            B = getattr(calc.values(self,i_self,t-tdel0,s,tele_angle_l=tele_l0,tele_angle_r=tele_angle,beam_angle_l=beam_angle_l,beam_angle_r=beam_angle_r,offset_l=offset_l0+offset_l,offset_r=offset_r0,mode='rec',ret=[ret2]),ret2)

        return B

    def twoPAAM_tele_aim_SS(self):
        # For Step-and_Stare
        tele_l_tot=[0,0,0]
        tele_r_tot=[0,0,0]
        tele_l_adjust=[0,0,0]
        tele_r_adjust=[0,0,0]
        t_l_adjust=[0,0,0]
        t_r_adjust=[0,0,0]
        
        if self.aimset.testSS==False:
            i_start = 1
            i_end = 4
        else:
            i_start=self.aimset.testSS
            i_end = i_start+1

        print('Only calculating for one laserlink')
        
        for link in range(i_start,i_end):
            t_plot = self.data.t_all[2:-3]
            t_adjust,[tele_l,tele_r],i_left,i_right = calc.SS_value(self,link,t_plot[0],t_plot[-1],'solve',self.aimset.FOV,ret=None,print_on=False,value=0)
            f_l = lambda t: calc.get_SS_func(t_adjust,tele_l,t)
            f_r = lambda t: calc.get_SS_func(t_adjust,tele_r,t)
            tele_l_tot[i_left-1] = f_l
            tele_r_tot[i_right-1] = f_r
            t_l_adjust[i_left-1] = t_adjust
            t_r_adjust[i_right-1] = t_adjust
            tele_l_adjust[i_left-1] = tele_l
            tele_r_adjust[i_right-1] = tele_r
        
        tele_l_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'l',x=t_l_adjust[i-1],y=tele_l_adjust[i-1])
        tele_r_ang = lambda i,t: calc.get_tele_SS(False,False,i,t,'r',x=t_r_adjust[i-1],y=tele_r_adjust[i-1])

        self.t_l_adjust = lambda i,t: t_l_adjust[i-1]
        self.t_r_adjust = lambda i,t: t_r_adjust[i-1]
        self.tele_l_ang_SS = lambda i,t: tele_l_ang
        self.tele_r_ang_SS = lambda i,t: tele_r_ang

        self.tele_l_ang_func = self.tele_l_ang_SS
        self.tele_r_ang_func = self.tele_r_ang_SS

        self.t_adjust = [t_l_adjust,t_r_adjust]
        self.tele_adjust = [tele_l_adjust,tele_r_adjust]
        self.tele_adjust_samp=[tele_l_adjust,tele_r_adjust]
        delat=['tele_l_ang','tele_r_ang','tele_l_ang_func','tele_r_ang_func','offset']
        for d in delat:
            try:
                delattr(self,d)
            except AttributeError:
                pass
        self.tele_l_ang = self.tele_l_ang_SS
        self.tele_r_ang = self.tele_r_ang_SS

        offset_calc={}
        offset_calc['l'] = lambda i,t: self.twoPAAM_PAAMin_aim_SS(i,t,'l')
        offset_calc['r'] = lambda i,t: self.twoPAAM_PAAMin_aim_SS(i,t,'r')
        self.offset = lambda i,t,s: offset_calc[s](i,t)

        return 0

    def twoPAAM_PAAMin_aim_SS(self,i,t,s):
        if s=='l':
            tdel = self.data.L_sl(i,t)
        elif s=='r':
            tdel = self.data.L_sr(i,t)

        if self.data.stat.calc_method=='Waluschka':
            tdel0 = tdel
        elif self.data.stat.calc_method=='Abram':
            tdel0 = 0.0

        ret = 'off'
        A = getattr(calc.values(self,i,t+tdel0,s,tele_angle_l=False,tele_angle_r=False,beam_angle_l=False,beam_angle_r=False,offset_l=0.0,offset_r=0.0,mode='send',ret=[ret]),ret)
        offset_calc = np.sign(A[2])*abs(np.arctan(A[2]/A[0]))

#        if s=='l':
#            print(getattr(calc.values(self,i,t+tdel0,s,tele_angle_l=False,tele_angle_r=False,beam_angle_l=False,beam_angle_r=False,offset_l=offset_calc,offset_r=0.0,mode='send',ret=[ret]),ret))

        return offset_calc


    def twoPAAM_tele_aim(self,i,t,s,test=False,**kwargs):
        #[self.tele_l_ang_fc,self.tele_r_ang_fc] = self.aimset.tele_control_ang_fc(option='center',value=0)
        
        if self.aimset.tele_control=='full_control' or test==True:
            ret = 'off'
            tele_l0 = -np.radians(30.0)
            tele_r0 = np.radians(30.0)
            
            [i_self,i_left,i_right] = const.i_slr(i)

            if s=='l':
                s_inv = 'r'
                tdel = self.data.L_rl(i,t)
                if self.data.stat.calc_method=='Waluschka':
                    tdel0 = tdel
                elif self.data.stat.calc_method=='Abram':
                    tdel0 = 0.0
                Al1 = getattr(calc.values(self,i_self,t+tdel0,s,tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=False,beam_angle_r=False,offset_l=0.0,offset_r=0.0,mode='send',ret=[ret]),ret)
                Ar1 = getattr(calc.values(self,i_left,t-tdel,s_inv,tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=False,beam_angle_r=False,offset_l=0.0,offset_r=0.0,mode='send',ret=[ret]),ret)
                offset_l1 = np.sign(Al1[2])*abs(np.arctan(Al1[2]/Al1[0]))
                offset_r1 = np.sign(Ar1[2])*abs(np.arctan(Ar1[2]/Ar1[0]))

            if s=='r':
                s_inv = 'l'
                tdel = self.data.L_rr(i,t)
                if self.data.stat.calc_method=='Waluschka':
                    tdel0 = tdel
                elif self.data.stat.calc_method=='Abram':          
                    tdel0 = 0.0
                Al1 = getattr(calc.values(self,i_right,t-tdel,s_inv,tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=False,beam_angle_r=False,offset_l=0.0,offset_r=0.0,mode='send',ret=[ret]),ret)
                Ar1 = getattr(calc.values(self,i_self,t+tdel0,s,tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=False,beam_angle_r=False,offset_l=0.0,offset_r=0.0,mode='send',ret=[ret]),ret)
                offset_l1 = np.sign(Al1[2])*abs(np.arctan(Al1[2]/Al1[0]))
                offset_r1 = np.sign(Ar1[2])*abs(np.arctan(Ar1[2]/Ar1[0]))

            ret2 = 'angx_wf_rec'
            
            if s=='l':
                tdel = self.data.L_rl(i,t)
                if self.data.stat.calc_method=='Waluschka':
                    tdel0 = tdel
                elif self.data.stat.calc_method=='Abram':
                    tdel0 = 0.0

                ang_extra_solve = -getattr(calc.values(self,i_self,t-tdel0,s,tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=False,beam_angle_r=False,offset_l=offset_l1,offset_r=offset_r1,mode='rec',ret=[ret2]),ret2)
            
            elif s=='r':
                tdel = self.data.L_rr(i,t)
                if self.data.stat.calc_method=='Waluschka':
                    tdel0 = tdel
                elif self.data.stat.calc_method=='Abram':
                    tdel0 = 0.0
                ang_extra_solve = - getattr(calc.values(self,i_self,t-tdel0,s,tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=False,beam_angle_r=False,offset_l=offset_l1,offset_r=offset_r1,mode='rec',ret=[ret2]),ret2)
                                
            if s=='l':
                ret = tele_l0+ang_extra_solve
                offset = offset_l1-ang_extra_solve

            elif s=='r':
                ret = tele_r0+ang_extra_solve
                offset = offset_r1-ang_extra_solve

        if self.aimset.tele_control=='no_control':
            if s=='l':
                ret = -np.radians(30.0)
                offset = 0.0

            elif s=='r':
                ret = np.radians(30.0)
                offset = 0.0

        return [ret,offset]
   
    def twoPAAM_PAAMout_aim(self,**kwargs):
        if self.aimset.PAAMout_control=='full_control':
            tele_l_ang0 = lambda i,t: -np.radians(30.0)
            tele_r_ang0 = lambda i,t: np.radians(30.0)
            [self.beam_l_ang_fc,self.beam_r_ang_fc] = self.PAAM_control_ang_fc(option='center',tele_l_ang = tele_l_ang0,tele_r_ang=tele_r_ang0)
        elif self.aimset.PAAMout_control=='no_control':
            self.beam_l_ang_fc = lambda i,t: 0.0
            self.beam_r_ang_fc = lambda i,t: 0.0

        return 0

 
