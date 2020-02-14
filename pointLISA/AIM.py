from imports import * 
import output
import numpy as np
import os
import yaml
from pointLISA import *
import pointLISA
# This class obtaines the pointing of the PAAM and telescope


class AIM():
    def __init__(self,data=False,setting=utils.Object(),filename=False,**kwargs):        
        from pointLISA import *
        # Get settings and parameters
        for key, value in kwargs.items():
            setattr(self,key,value)

        if data!=None:
            self.data = data
            aimset0 = setting
            aimset = utils.Object()
                 
            for key,value in kwargs.items():
                if filename!=False:
                    try:
                        getattr(aimset,key)
                    except:
                        setattr(self,key,value)
                        setattr(aimset,key,value)
                else:
                        setattr(self,key,value)
                        setattr(aimset,key,value)
            
            for key in setting.__dict__.keys():
                if key not in aimset.__dict__.keys():
                    setattr(aimset,key,setting.__dict__[key])
                    setattr(self,key,setting.__dict__[key])

            for k in settings.aimset.__dict__.keys():
                if k not in aimset.__dict__.keys():
                    setattr(aimset,k,settings.aimset.__dict__[k])
            self.aimset = aimset

            

            if data!=False:
                if self.data.input_file==None:
                    print('Start calculating telescope and PAAM aim')
                    self.get_offset_inplane(self.aimset.offset_tele)
                

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
                ret = methods.read(direct=option)
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
                    read_folder = os.path.dirname(os.path.realpath(__file__))+'/parameters/'+self.data.calc_method+'/'

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
                            pos_l = methods.get_nearest_smaller_value(keys,days)
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


    def tele_control_ang_fc(self,option=None,value=False):
        '''Obtains the telescope pointing angles for a continuous actuation (full_control)'''
        # Option 'wavefront' means poiting with the purpose of getting a zero/small tilt of the receiving wavefront
        # 'center' means pointing it to te center of the receiving telescope aperture
        if option==None:
            option = self.aimset.tele_option

        print('Telescope pointing strategy: '+option)
        
        max_count=5
        scale=1
        
        ang_l = lambda i,t: methods.tele_point_calc(self,i,t,'l',option,max_count=max_count,scale=scale,value=value,tele_l0=np.radians(-30.0),tele_r0=np.radians(30.0),beam_l0=0.0,beam_r0=0.0)
        ang_r = lambda i,t: methods.tele_point_calc(self,i,t,'r',option,max_count=max_count,scale=scale,value=value,tele_l0=np.radians(-30.0),tele_r0=np.radians(30.0),beam_l0=0.0,beam_r0=0.0)
        
        self.option_tele = option
        self.tele_control = 'full_control'

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

        if self.data.input_file==None:
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
                
                if option=='wavefront':
                    ret = 'angx_wf_send'
                    lim=self.aimset.FOV #Klopt niet want alleen in inplane i.p.v. totaal, kan ook met I ...adjust
                elif option=='center':
                    ret = 'xoff'
                    lim=self.aimset.width/2.0

                for link in range(1,4):
                    t_plot = self.data.t_all[2:-3]
                    t_adjust,[tele_l,tele_r],i_left,i_right = methods.SS_value(self.aim0,link,t_plot[0],t_plot[-1],'solve',lim,ret=ret,print_on=False,value=value)
                    f_l = lambda t: methods.get_SS_func(t_adjust,tele_l,t)
                    f_r = lambda t: methods.get_SS_func(t_adjust,tele_r,t)
                    tele_l_tot[i_left-1] = f_l
                    tele_r_tot[i_right-1] = f_r
                    t_l_adjust[i_left-1] = t_adjust
                    t_r_adjust[i_right-1] = t_adjust

                self.t_l_adjust = lambda i,t: t_l_adjust[i-1]
                self.t_r_adjust = lambda i,t: t_r_adjust[i-1]
                self.tele_l_ang_SS = lambda i,t: tele_l_tot[i-1](t)
                self.tele_r_ang_SS = lambda i,t: tele_r_tot[i-1](t)

                self.tele_l_ang_func = self.tele_l_ang_SS
                self.tele_r_ang_func = self.tele_r_ang_SS

                self.t_adjust = [t_l_adjust,t_r_adjust]
                self.tele_adjust = [tele_l_tot,tele_r_tot]
                self.tele_adjust_samp=[tele_l,tele_r]

            elif type(method)==list and method[0]=='Imported pointing':
                # Import (previously obtained) angles
                print(method[0])
                self.tele_l_ang = lambda i,t: methods.get_tele_SS(False,False,i,t,'l',x=method[1]['SC'+str(i)+', left']['x'],y=method[1]['SC'+str(i)+', left']['y'])
                self.tele_r_ang = lambda i,t: methods.get_tele_SS(False,False,i,t,'r',x=method[1]['SC'+str(i)+', right']['x'],y=method[1]['SC'+str(i)+', right']['y'])

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

            else:
                raise ValueError('Please select valid telescope pointing method')
    
        else:
            print('Importing pointing angles from:')
            print(self.data.input_file)
            ret = pointLISA.read_write.read_output(filenames=self.data.input_file)
            tele_l_ang=[]
            tele_r_ang=[]
            if method=='SS':
                self.tele_l_ang = lambda i,t: methods.get_tele_SS(False,False,i,t,'l',x=ret['SC'+str(i)+', left']['x'],y=ret['SC'+str(i)+', left']['y'])
                self.tele_r_ang = lambda i,t: methods.get_tele_SS(False,False,i,t,'r',x=ret['SC'+str(i)+', right']['x'],y=ret['SC'+str(i)+', right']['y'])

                t_adjust={}
                tele_ang_adjust = {}
                for i in range(1,4):
                    t_adjust[str(i)]={}
                    tele_ang_adjust[str(i)]={}
                    
                    [t_adjust[str(i)]['l'],tele_ang_adjust[str(i)]['l']]=getattr(ret[0].l,'i'+str(i)).adjust[0]
                    [t_adjust[str(i)]['r'],tele_ang_adjust[str(i)]['r']]=getattr(ret[0].r,'i'+str(i)).adjust[0]

                self.tele_l_ang = lambda i,t: methods.get_tele_SS(False,False,i,t,'l',x=t_adjust[str(i)]['l'],y=tele_ang_adjust[str(i)]['l'])
                self.tele_r_ang = lambda i,t: methods.get_tele_SS(False,False,i,t,'l',x=t_adjust[str(i)]['r'],y=tele_ang_adjust[str(i)]['r'])
                
                self.t_adjust = t_adjust
                self.tele_ang_adjust = tele_ang_adjust

            else:

                for i in range(1,4):
                    t_l =  getattr(getattr(ret[0],'l'),'i'+str(i)).tele_ang
                    t_r =  getattr(getattr(ret[0],'r'),'i'+str(i)).tele_ang
                    tele_l_ang.append(methods.interpolate(t_l[0][0],t_l[0][1]))
                    tele_r_ang.append(methods.interpolate(t_r[0][0],t_r[0][1]))

                self.tele_l_ang = lambda i,t: tele_l_ang[i-1](t)
                self.tele_r_ang = lambda i,t: tele_r_ang[i-1](t)
            
            if self.PAAM_deg==2:
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
                    offset['l'][i] = methods.interpolate(offset_l[0][0],offset_l[0][1])
                    offset['r'][i] = methods.interpolate(offset_r[0][0],offset_r[0][1])
                self.offset={}
                self.offset = offset

            else:
                self.get_offset_inplane(self.aimset.offset_tele)

        try:
            self.tele_l_ang
        except AttributeError:
            if self.sampled==True:
                try:
                    self.t_sample
                except AttributeError:
                    self.t_sample = methods.get_t_sample(self,speed=self.aimset.sample_speed)
                
                tele_l_ang=[]
                tele_r_ang=[]
                print("Sampling and fitting telescope angles") #...add offset voor PAAM_deg==2
                for i in range(1,4):
                    tele_l_ang.append(methods.interpolate(self.t_sample['l'][i-1],np.array([self.tele_l_ang_func(i,t) for t in self.t_sample['l'][i-1]])))
                    tele_r_ang.append(methods.interpolate(self.t_sample['r'][i-1],np.array([self.tele_r_ang_func(i,t) for t in self.t_sample['r'][i-1]])))
                self.tele_l_ang_samp = lambda i,t: tele_l_ang[i-1](t)
                self.tele_r_ang_samp = lambda i,t: tele_r_ang[i-1](t)
                
                self.tele_l_ang = self.tele_l_ang_samp
                self.tele_r_ang = self.tele_r_ang_samp

            else:
                self.tele_l_ang = self.tele_l_ang_func
                self.tele_r_ang = self.tele_r_ang_func

        return 0

    def PAAM_control_ang_fc(self,option='center'):
        '''Obtains the PAAM pointing angles for a continuous actuation (full_control)'''
        print('PAAM pointing strategy: '+option)

        if option=='center':
            ang_l = lambda i,t: output.PAAM_center_calc(self,i,t,tele_l=self.tele_l_ang(i,t),tele_r=self.tele_r_ang(utils.i_slr(i)[1],t),lim=self.aimset.limit_yoff,method=self.aimset.PAAM_method_solve,para=self.aimset.optimize_PAAM,value=self.aimset.optimize_PAAM_value,margin=self.aimset.optimize_PAAM_margin)[0][0]
            ang_r = lambda i,t: output.PAAM_center_calc(self,utils.i_slr(i)[2],t,tele_l=self.tele_l_ang(utils.i_slr(i)[2],t),tele_r=self.tele_r_ang(i,t),lim=self.aimset.limit_yoff,method=self.aimset.PAAM_method_solve,para=self.aimset.optimize_PAAM,value=self.aimset.optimize_PAAM_value,margin=self.aimset.optimize_PAAM_margin)[0][1]

        elif option=='wavefront':
            ang_l = lambda i,t: output.PAAM_wavefront_calc(self,i,t,'l',lim=self.aimset.limit_angy)
            ang_r = lambda i,t: output.PAAM_wavefront_calc(self,i,t,'r',lim=self.aimset.limit_angy)

        self.PAAM_control = 'full_control'
        self.option_PAAM = option
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
        
        if self.data.input_file==None:
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
                
                if self.sampled==True:
                    sampled = self.sample()
                    self.aim_sampled = sampled
                else:
                    self.aim_sampled=False

            if self.inp!=False:
                self.aim_sampled = self

        else:
            print('Importing pointing angles from:')
            print(self.data.input_file)
            ret = pointLISA.read_write.read_output(filenames=self.data.input_file)
            beam_l_ang=[]
            beam_r_ang=[]
            for i in range(1,4):
                b_l =  getattr(getattr(ret[0],'l'),'i'+str(i)).PAAM_ang
                b_r =  getattr(getattr(ret[0],'r'),'i'+str(i)).PAAM_ang
                beam_l_ang.append(methods.interpolate(b_l[0][0],b_l[0][1]))
                beam_r_ang.append(methods.interpolate(b_r[0][0],b_r[0][1]))

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
            self.t_sample = methods.get_t_sample(self,speed=self.aimset.sample_speed)

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
                tele_l_ang.append(methods.interpolate(self.t_sample['l'][i-1],np.array([self.tele_l_ang(i,t) for t in self.t_sample['l'][i-1]])))
                tele_r_ang.append(methods.interpolate(self.t_sample['r'][i-1],np.array([self.tele_r_ang(i,t) for t in self.t_sample['l'][i-1]])))
            beam_l_ang.append(methods.interpolate(self.t_sample['l'][i-1],np.array([self.beam_l_ang(i,t) for t in self.t_sample['l'][i-1]])))
            beam_r_ang.append(methods.interpolate(self.t_sample['r'][i-1],np.array([self.beam_r_ang(i,t) for t in self.t_sample['r'][i-1]])))
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
        tele_l_coor = methods.coor_tele(self.data,i,t,tele_l_ang(i,t))
        tele_r_coor = methods.coor_tele(self.data,i,t,tele_r_ang(i,t))
        tele_l_vec = LA.unit(tele_l_coor[0])*self.data.L_tele
        tele_r_vec = LA.unit(tele_r_coor[0])*self.data.L_tele
        tele_l_start = tele_l_vec+np.array(self.data.putp(i,t))
        tele_r_start = tele_r_vec+np.array(self.data.putp(i,t))

        return [[tele_l_coor,tele_r_coor],[tele_l_vec,tele_r_vec],[tele_l_start,tele_r_start]]

    def get_beam_coor(self,i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang,offset): #...only send frame
        ''' Calculating new pointing vectors and beam coordinate system '''
        if offset == False:
            offset = self.offset_tele

        offset_l = methods.get_offset(self,i,t,'l')
        offset_r = methods.get_offset(self,i,t,'r')

        beam_l_coor = methods.beam_coor_out(self.data,i,t,tele_l_ang(i,t),beam_l_ang(i,t),offset_l)
        beam_r_coor = methods.beam_coor_out(self.data,i,t,tele_r_ang(i,t),beam_r_ang(i,t),offset_r)

        # Calculating the Transmitted beam direction and position of the telescope aperture
        beam_l_direction = beam_l_coor[0]
        beam_r_direction = beam_r_coor[0]
        beam_l_start = self.data.L_tele*beam_l_direction+np.array(self.data.putp(i,t)) #...kan weg
        beam_r_start = self.data.L_tele*beam_r_direction+np.array(self.data.putp(i,t)) #...kan weg

        return [[beam_l_coor,beam_r_coor],[beam_l_direction,beam_r_direction],[beam_l_start,beam_r_start]]

    def get_coordinate_systems(self,iteration_val=False,option='self'):
        '''Make functions of the coordinate vectors and systems'''
        if iteration_val==False:
            tele_l_ang = self.tele_l_ang
            tele_r_ang = self.tele_r_ang
            beam_l_ang = self.beam_l_ang
            beam_r_ang = self.beam_r_ang
            offset = self.offset_tele
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

        [i_self,i_left,i_right] = utils.i_slr(i)
        if side=='l':
            if self.tele_control=='full_control':
                Dt = self.data.L_sl_func_tot(i_self,t)
                ang_in_l = methods.tele_point_calc(self,i,t,'l','center',max_count=max_count,scale=scale,value=value,tele_l0=np.radians(-30.0),tele_r0=np.radians(30.0),beam_l0=0.0,beam_r0=0.0,offset_l0=0.0,offset_r0=0.0)
                ang_in_r = methods.tele_point_calc(self,utils.i_slr(i)[1],t+Dt,'r','center',max_count=max_count,scale=scale,value=value,tele_l0=np.radians(-30.0),tele_r0=np.radians(30.0),beam_l0=0.0,beam_r0=0.0,offset_l0=0.0,offset_r0=0.0)
            elif self.tele_control=='no_control':
                ang_in_l = np.radians(-30.0)
                ang_in_r = np.radians(30.0)

            if self.PAAM_control=='full_control':
                ang_out = output.PAAM_center_calc(self,i,t,tele_l=ang_in_l,tele_r=ang_in_r,lim=self.aimset.limit_yoff,method=self.aimset.PAAM_method_solve,para=self.aimset.optimize_PAAM,value=self.aimset.optimize_PAAM_value,margin=self.aimset.optimize_PAAM_margin,beam_l=0.0,beam_r=0.0,offset_l=0.0,offset_r=0.0)[0][0]
            elif self.PAAM_control=='no_control':
                ang_out = 0.0

            return [ang_in_l,ang_out]

        elif side=='r':
            if self.tele_control=='full_control':
                Dt = self.data.L_sr_func_tot(i_self,t)
                ang_in_l = methods.tele_point_calc(self,utils.i_slr(i)[2],t+Dt,'l','center',max_count=max_count,scale=scale,value=value,tele_l0=np.radians(-30.0),tele_r0=np.radians(30.0),beam_l0=0.0,beam_r0=0.0,offset_l0=0.0,offset_r0=0.0)
                ang_in_r = methods.tele_point_calc(self,i,t,'r','center',max_count=max_count,scale=scale,value=value,tele_l0=np.radians(-30.0),tele_r0=np.radians(30.0),beam_l0=0.0,beam_r0=0.0,offset_l0=0.0,offset_r0=0.0)
            elif self.tele_control=='no_control':
                ang_in_l=np.radians(-30.0)
                ang_in_r=np.radians(30.0)

            if self.PAAM_control=='full_control':
                ang_out = output.PAAM_center_calc(self,utils.i_slr(i)[2],t,tele_l=ang_in_l,tele_r=ang_in_r,lim=self.aimset.limit_yoff,method=self.aimset.PAAM_method_solve,para=self.aimset.optimize_PAAM,value=self.aimset.optimize_PAAM_value,margin=self.aimset.optimize_PAAM_margin,beam_l=0.0,beam_r=0.0,offset_l=0.0,offset_r=0.0)[0][1]
            elif self.PAAM_control=='no_control':
                ang_out = 0.0

            return [ang_in_r,ang_out]
    
    def twoPAAM_pointing(self,i,t,side,out,mode,short=False):
        '''Obtains the dual axis PAAM poiting angles (which uses the twoPAAM_calc helper function'''
        [i_self,i_left,i_right] = utils.i_slr(i)
        if mode=='rec':
            if side=='l':
                Dt = self.data.L_rl_func_tot(i_self,t)
                [ang_in_start,ang_out_start] = self.twoPAAM_calc(i_left,t-Dt,'r')
                [ang_in_end,ang_out_end] = self.twoPAAM_calc(i_self,t,'l')
                i_rec = i_self
                i_send = i_left
            elif side=='r':
                Dt = self.data.L_rr_func_tot(i_self,t)
                [ang_in_start,ang_out_start] = self.twoPAAM_calc(i_right,t-Dt,'l')
                [ang_in_end,ang_out_end] = self.twoPAAM_calc(i_self,t,'r')
                i_rec = i_self
                i_send = i_right
            if short==False and self.offset_tele!='no_control':
                coor_tele_rec = methods.coor_tele(self.data,i_rec,t,ang_in_end)
                coor_tele_send = methods.coor_tele(self.data,i_send,t-Dt,ang_in_start)
                pos_tele_rec = LA.unit(coor_tele_rec[0])*self.data.L_tele+np.array(self.data.putp(i_rec,t))
                coor_beam_send = methods.beam_coor_out(self.data,i_send,t-Dt,ang_in_start,ang_out_start,0.0)
                pos_tele_send = LA.unit(coor_tele_send[0])*self.data.L_tele+np.array(self.data.putp(i_send,t-Dt))

        elif mode=='send':
            if side=='l':
                Dt = self.data.L_sl_func_tot(i_self,t)
                [ang_in_start,ang_out_start] = self.twoPAAM_calc(i_self,t,'l')
                [ang_in_end,ang_out_end] = self.twoPAAM_calc(i_left,t+Dt,'r')
                i_rec = i_left
                i_send = i_self
            elif side=='r':
                Dt = self.data.L_sr_func_tot(i_self,t)
                [ang_in_start,ang_out_start] = self.twoPAAM_calc(i_self,t,'r')
                [ang_in_end,ang_out_end] = self.twoPAAM_calc(i_right,t+Dt,'l')
                i_rec = i_right
                i_send = i_self
            if short==False and self.offset_tele!='no_control':
                coor_tele_rec = methods.coor_tele(self.data,i_rec,t+Dt,ang_in_end)
                coor_tele_send = methods.coor_tele(self.data,i_send,t,ang_in_start)
                pos_tele_rec = LA.unit(coor_tele_rec[0])*self.data.L_tele+np.array(self.data.putp(i_rec,t+Dt))
                coor_beam_send = methods.beam_coor_out(self.data,i_send,t,ang_in_start,ang_out_start,0.0)
                pos_tele_send = LA.unit(coor_tele_send[0])*self.data.L_tele+np.array(self.data.putp(i_send,t))
        
        if short==False and self.offset_tele!='no_control':
            z = np.linalg.norm(LA.matmul(coor_beam_send,pos_tele_rec - pos_tele_send))

            R = out.R(z)
            R_vec = -R*coor_beam_send[0]
            R_vec_tele_rec = LA.unit(LA.matmul(coor_tele_rec,-R_vec))
            offset =  np.sign(R_vec_tele_rec[2])*np.arctan(abs(R_vec_tele_rec[2]/R_vec_tele_rec[0]))
            [z_off,y_off,x_off] = LA.matmul(coor_beam_send,pos_tele_rec - pos_tele_send) #are correct

        elif short==True:
            offset=np.nan
        elif self.offset_tele=='no_control':
            offset=0.0
           
        return [ang_in_end,ang_out_end,ang_in_start,ang_out_start,t,Dt,i_rec,i_send,offset]

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

        self.twoPAAM_tele_aim()
        self.twoPAAM_PAAMout_aim()
        self.twoPAAM_PAAMin_aim()

        if (self.sampled == True and sampled==None) or sampled==True:
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
                
                tele_func_dict['l'][i] = methods.interpolate(np.array(t_l_calc),np.array(tele_l_calc))
                beam_func_dict['l'][i] = methods.interpolate(np.array(t_l_calc),np.array(beam_l_calc))
                tele_func_dict['r'][i] = methods.interpolate(np.array(t_r_calc),np.array(tele_r_calc))
                beam_func_dict['r'][i] = methods.interpolate(np.array(t_r_calc),np.array(beam_r_calc))
#
                offset['l'][i] = methods.interpolate(np.array(t_l_calc),np.array(offset_l_calc))
                offset['r'][i] = methods.interpolate(np.array(t_r_calc),np.array(offset_r_calc))
            
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


    def twoPAAM_get_tele(self,i,t,side,out,short=False):
        '''Getfunction to obtain the telescope pointing angle for a dual axis PAAM pointing'''
        A = self.twoPAAM_pointing(i,t,side,out,'rec',short=short)
        return A[0] - A[-1]
    
    def twoPAAM_get_beam(self,i,t,side,out,short=False):
        '''Getfunction to obtain the beam direction for a dual axis PAAM pointing'''
        A = self.twoPAAM_pointing(i,t,side,out,'rec',short=short)
        return A[1]

    def twoPAAM_get_offset(self,i,t,side,out,short=False):
        '''Getfunction to obtain the telescopei offset angle for a dual axis PAAM pointing'''
        A = self.twoPAAM_pointing(i,t,side,out,'rec',short=short)
        return A[-1]

    def twoPAAM_tele_aim(self,**kwargs):
        '''Obtianing the telescope pointing angle when using a dual axis PAAM'''
        try:
            tele_control = kwargs['tele_control']
            delattr(self,'tele_control')
            setattr(self,'tele_control',tele_control)
        except KeyError:
            pass
        print('tele_control: '+self.tele_control)

        if self.tele_control=='no_control':
            tele_l_ang = lambda i,t: np.radians(-30.0)
            tele_r_ang = lambda i,t: np.radians(30.0)
        
        elif self.tele_control=='full_control':
            ret='angx_arm_tele_rec'
            tele_l0 = np.radians(-30.0)
            tele_r0 = np.radians(30.0)
            tele_l_ang = lambda i,t: tele_l0+getattr(output.values(self,i,t,'l',tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=0.0,beam_angle_r=0.0,offset_l=0.0,offset_r=0.0,mode='rec',ret=[ret]),ret)
            tele_r_ang = lambda i,t: tele_r0+getattr(output.values(self,i,t,'r',tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=0.0,beam_angle_r=0.0,offset_l=0.0,offset_r=0.0,mode='rec',ret=[ret]),ret)
        
        elif self.tele_control=='offset':
            ret='angx_arm_tele_rec'
            tele_l0 = np.radians(-30.0)
            tele_r0 = np.radians(30.0)
            try:
                self.t_sample
            except AttributeError:
                self.t_sample = methods.get_t_sample(self,speed=self.aimset.sample_speed)
            
            tele_l_samp={}
            tele_r_samp={}
            for i in range(1,4):
                tele_l_samp[i]=[]
                tele_r_samp[i]=[]
                for t in self.t_sample['l'][1]:
                    tele_l_samp[i].append(tele_l0+getattr(output.values(self,i,t,'l',tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=0.0,beam_angle_r=0.0,offset_l=0.0,offset_r=0.0,mode='rec',ret=[ret]),ret))
                    tele_r_samp[i].append(tele_r0+getattr(output.values(self,i,t,'r',tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=0.0,beam_angle_r=0.0,offset_l=0.0,offset_r=0.0,mode='rec',ret=[ret]),ret))
            
            tele_l_ang = lambda i,t: np.nanmean(tele_l_samp[i])
            tele_r_ang = lambda i,t: np.nanmean(tele_r_samp[i])

        elif self.tele_control=='SS':
            try:
                t0 = kwargs['t0']
            except KeyError:
                t0 = self.data.t_all[1]
            try:
                t_end = kwargs['t_end']
            except KeyError:
                t_end = self.data.t_all[-1]
            method='solve'
            lim=self.aimset.FOV
            print('lim= ',lim)
            dt=24*3600
            scale = self.aimset.tele_SS_scale
            SS=[]

            for link in range(1,4):
                print('Link is '+str(link))
                SS.append(methods.SS_value(self,link,t0,t_end,method,lim,print_on=True,dt=dt,scale=scale))

            tele_l_adjust={}
            tele_r_adjust={}
            t_l_adjust={}
            t_r_adjust={}
            for j in range(0,len(SS)):
                tele_l_adjust[SS[j][2]] = SS[j][1][0]
                tele_r_adjust[SS[j][3]] = SS[j][1][1]
                t_l_adjust[SS[j][2]] = SS[j][0]
                t_r_adjust[SS[j][3]] = SS[j][0]
            
            tele_l_ang = lambda i,t: methods.get_SS_func(t_l_adjust[i],tele_l_adjust[i],t)
            tele_r_ang = lambda i,t: methods.get_SS_func(t_r_adjust[i],tele_r_adjust[i],t)
            tele_adjust={'l': tele_l_adjust,'r':tele_r_adjust}
            t_adjust={'l': t_l_adjust,'r':t_r_adjust}
        else: 
            raise ValueError('Please select a valid telescope pointing method (tele_control)')
        
        param=['tele_l_ang','tele_r_ang','t_l_adjust','t_r_adjust','tele_l_adjust','tele_r_adjust','tele_adjust','t_adjust']
        for p in param:
            try:
                delattr(self,p)
            except AttributeError:
                pass
            try:
                setattr(self,p,locals()[p])
                print('New '+p+' value')
            except KeyError:
                pass

        return 0
        
    def twoPAAM_PAAMout_aim(self,**kwargs):
        '''Runner function to get the outpane PAAM poiting angle of the dual axis PAAM'''
        try:
            PAAMiout_control = kwargs['PAAMout_control']
            delattr(self,'PAAMout_control')
            setattr(self,'PAAMout_control',PAAMout_control)
        except KeyError:
            pass
        print('PAAMout_control: '+self.PAAMout_control)

        if self.PAAMout_control=='no_control':
            beam_l_ang = lambda i,t: 0.0
            beam_r_ang = lambda i,t: 0.0

        elif self.PAAMout_control=='full_control':
            tele_l0=False
            tele_r0=False
            ret = 'angy_arm_tele_send'

            beam_l_ang = lambda i,t: -getattr(output.values(self,i,t,'l',tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=0.0,beam_angle_r=0.0,offset_l=0.0,offset_r=0.0,mode='send',ret=[ret]),ret)
            beam_r_ang = lambda i,t: -getattr(output.values(self,i,t,'r',tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=0.0,beam_angle_r=0.0,offset_l=0.0,offset_r=0.0,mode='send',ret=[ret]),ret)
        
        elif self.PAAMout_control=='offset':
            tele_l0=False
            tele_r0=False
            ret = 'angy_arm_tele_send'
            try:
                self.t_sample
            except AttributeError:
                self.t_sample = methods.get_t_sample(self,speed=self.aimset.sample_speed)

            beam_l_samp={}
            beam_r_samp={}
            for i in range(1,4):
                beam_l_samp[i]=[]
                beam_r_samp[i]=[]
                for t in self.t_sample['l'][1]:
                    beam_l_samp[i].append(-getattr(output.values(self,i,t,'l',tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=0.0,beam_angle_r=0.0,offset_l=0.0,offset_r=0.0,mode='send',ret=[ret]),ret))
                    beam_r_samp[i].append(-getattr(output.values(self,i,t,'r',tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=0.0,beam_angle_r=0.0,offset_l=0.0,offset_r=0.0,mode='send',ret=[ret]),ret))
            
            beam_l_ang = lambda i,t: np.nanmean(beam_l_samp[i])
            beam_r_ang = lambda i,t: np.nanmean(beam_r_samp[i])

        else:
            raise ValueError('Please select a valid PAAM outplane pointing method (PAAMout_control)')

        param=['beam_l_ang','beam_r_ang']
        for p in param:
            try:
                delattr(self,p)
            except AttributeError:
                pass
            try:
                setattr(self,p,locals()[p])
                print('New '+p+' value')
            except KeyError:
                pass

        return 0

    def twoPAAM_PAAMin_aim(self,**kwargs):
        '''Runner function to get the inplane PAAM poiting angle of th
e dual axis PAAM'''
        try:
            PAAMin_control = kwargs['PAAMin_control']
            delattr(self,'PAAMin_control')
            setattr(self,'PAAMin_control',PAAMin_control)
        except KeyError:
            pass
        print('PAAMin_control: '+self.PAAMin_control)

        if self.PAAMin_control=='no_control':
            offset_calc={}
            offset_calc['l'] = lambda i,t: 0.0
            offset_calc['r'] = lambda i,t: 0.0
            offset=lambda i,t,s: offset_calc[s](i,t)

        elif self.PAAMin_control=='full_control':
            tele_l0=False
            tele_r0=False
            beam_l0=0.0
            beam_r0=0.0
            offset_calc={}
            ret = 'angx_arm_tele_send'

            offset_calc['l'] = lambda i,t: getattr(output.values(self,i,t,'l',tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=beam_l0,beam_angle_r=beam_r0,offset_l=0.0,offset_r=0.0,mode='send',ret=[ret]),ret)
            offset_calc['r'] = lambda i,t: getattr(output.values(self,i,t,'r',tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=beam_l0,beam_angle_r=beam_r0,offset_l=0.0,offset_r=0.0,mode='send',ret=[ret]),ret)
            
            offset=lambda i,t,s: offset_calc[s](i,t)
        
        elif self.PAAMin_control=='offset':
            tele_l0=False
            tele_r0=False
            beam_l0=0.0
            beam_r0=0.0
            ret = 'angx_arm_tele_send'
            offset_calc={}

            try:
                self.t_sample
            except AttributeError:
                self.t_sample = methods.get_t_sample(self,speed=self.aimset.sample_speed)
            
            offset_l_samp={}
            offset_r_samp={}
            for i in range(1,4):
                offset_l_samp[i]=[]
                offset_r_samp[i]=[]
                for t in self.t_sample['l'][1]:
                    offset_l_samp[i].append(getattr(output.values(self,i,t,'l',tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=beam_l0,beam_angle_r=beam_r0,offset_l=0.0,offset_r=0.0,mode='send',ret=[ret]),ret))
                    offset_r_samp[i].append(getattr(output.values(self,i,t,'r',tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_angle_l=beam_l0,beam_angle_r=beam_r0,offset_l=0.0,offset_r=0.0,mode='send',ret=[ret]),ret))


            offset_calc['l'] = lambda i,t: np.nanmean(offset_l_samp[i])
            offset_calc['r'] = lambda i,t: np.nanmean(offset_r_samp[i])
            offset=lambda i,t,s: offset_calc[s](i,t)

        else:
            raise ValueError('Please select a valid PAAM inplane pointing method (PAAMin_control)')

        param=['offset']
        for p in param:
            try:
                delattr(self,p)
            except AttributeError:
                pass
            try:
                setattr(self,p,locals()[p])
                print('New '+p+' value')
            except KeyError:
                pass

        return 0


