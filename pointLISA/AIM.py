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
        #import pointLISA
        from pointLISA import *
        for key, value in kwargs.items():
            setattr(self,key,value)

        if data!=None:
            self.data = data
            aimset0 = setting
            aimset = utils.Object()
            #self.aimset=settings
            #for key in settings.__dict__.keys():
            #    setattr(self,key,settings.__dict__[key])
            
            #if data.input_file!=False:
            #    aimset_read = read_write.read_options(data.input_file)
                 
            for key,value in kwargs.items():
                #print(key)
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
                print('Start calculating telescope and PAAM aim')
                self.get_offset_inplane(self.offset_tele)
                
                if self.init!=True:
                    try:
                        self.aim0 = kwargs['aim0']
                    except KeyError:
                        self.aim0=AIM.AIM(data=None)
                    self.angles0=kwargs.pop('angles0',False)
                    self.aim0.data=data
                    self.aim0.tele_l_ang = self.angles0[0]
                    self.aim0.tele_r_ang = self.angles0[2]
                    self.aim0.beam_l_ang = self.angles0[1]
                    self.aim0.beam_r_ang = self.angles0[3]
                    self.aim0.get_coordinate_systems(option='self')
                    self.aim0.offset_tele = self.offset_tele
                    for key in settings.__dict__.keys():
                        setattr(self.aim0,key,settings.__dict__[key])
                    self.aim0.get_offset_inplane(self.aim0.offset_tele)
                    setattr(self.aim0,'aimset',self.aimset)
                    

                    #self.angles_old=kwargs.pop('angles_old',False)
                    #self.aim_old=AIM(False)
                    #self.aim_old.data=data
                    #self.aim_old.tele_l_ang = self.angles_old[0]
                    #self.aim_old.tele_r_ang = self.angles_old[2]
                    #self.aim_old.beam_l_ang = self.angles_old[1]
                    #self.aim_old.beam_r_ang = self.angles_old[3]
                    #self.aim_old.offset_tele = self.offset_tele
                    #self.aim_old.get_coordinate_systems(option='self')

    def get_offset_inplane(self,option,**kwargs):
        offset = {'l': {1: 0.0, 2: 0.0, 3: 0.0},
 'r': {1: 0.0, 2: 0.0, 3: 0.0}}
        
        if option==0:
            pass

        elif option==True:
            offset = {'l': {1: -0.00018360462896226676, 2: 0.0, 3: 0.0},
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
                        out=out+line.split('\n')[0]
                offset = yaml.load(out)
                read_file.close()
        
        elif type(option==dict):
            offset = option

        else:
            raise ValueError("Please select offset tele values or method")

        self.offset_tele = offset

        return 0


    def tele_control_ang_fc(self,option='wavefront',value=False):
        # Option 'wavefront' means poiting with the purpose of getting a zero/small tilt of the receiving wavefront
        # 'center' means pointing it to te center of the receiving telescope aperture
        print('Telescope pointing strategy: '+option)
        
        max_count=5
        scale=1

        ang_l = lambda i,t: methods.tele_point_calc(self.aim0,i,t,'l',option,max_count=max_count,scale=scale,value=value)
        ang_r = lambda i,t: methods.tele_point_calc(self.aim0,i,t,'r',option,max_count=max_count,scale=scale,value=value)
        

        self.tele_option = option
        return [ang_l,ang_r]


    def tele_aim(self,method=False,lim=1e-10,tele_ang_extra=True,option=False):
        if method==False:
            method=self.tele_method
        else:
            self.tele_method=method
        
        if option=='wavefront':
            value=self.aimset.value_wavefront
        elif option=='center':
            value=self.aimset.value_center


        try:
            print('The telescope control method is: '+method)
        except:
            print('The telescope control method is: user defined')

        print(' ')

        #if self.inp==False:
        if self.data.input_file==None or self.init==True:
            if method=='no_control':
                if tele_ang_extra==False:
                    offset_l = {'1':0,'2':0,'3':0}
                    offset_r = {'1':0,'2':0,'3':0}
                elif tele_ang_extra==True:
                    try:
                        offset_l = self.wfe.mean_angin_l
                        offset_r = self.wfe.mean_angin_r
                    except:
                        self.wfe.do_mean_angin()
                        offset_l = self.wfe.mean_angin_l
                        offset_r = self.wfe.mean_angin_r
                else:
                    [offset_l,offset_r] = tele_ang_extra

                self.offset =[offset_l,offset_r]
                self.tele_l_ang_func = lambda i,t: np.radians(-30)+offset_l[str(i)]
                self.tele_r_ang_func = lambda i,t: np.radians(30)+offset_r[str(i)]

            elif method=='full_control':
                [self.tele_ang_l_fc,self.tele_ang_r_fc] = self.tele_control_ang_fc(option=option,value=value)
                self.tele_l_ang_func = self.tele_ang_l_fc
                self.tele_r_ang_func = self.tele_ang_r_fc

            elif 'min' in method:
                if 'spot'==method.split('_')[-1]:
                    self.tele_l_ang_func = lambda i,t: methods.spotsize_limit(self.data,self.aim0,i,t,'l',limit=-self.wfe.spotsize)
                    self.tele_r_ang_func = lambda i,t: methods.functions.spotsize_limit(self.data,self.aim0,i,t,'r',limit=-self.wfe.spotsize)
            elif 'max' in method:
                if 'spot'==method.split('_')[-1]:
                    self.tele_l_ang_func = lambda i,t: methods.spotsize_limit(self.data,self.aim0,i,t,'l',limit=self.wfe.spotsize)
                    self.tele_r_ang_func = lambda i,t: methods.functions.spotsize_limit(self.data,self.aim0,i,t,'r',limit=self.wfe.spotsize)

            elif 'SS' in method:
                tele_l_tot=[0,0,0]
                tele_r_tot=[0,0,0]
                t_l_adjust=[0,0,0]
                t_r_adjust=[0,0,0]

                for link in range(1,4):
                    t_plot = self.data.t_all[2:-3]
                    t_adjust,[tele_l,tele_r],i_left,i_right = methods.SS_value(self.aim0,link,t_plot[0],t_plot[-1],'solve',self.aimset.width/2.0,print_on=False,value=value)
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

            elif type(method)==list and method[0]=='Imported pointing':
                print(method[0])
                self.tele_l_ang = lambda i,t: pack.functions.get_tele_SS(False,False,i,t,'l',x=method[1]['SC'+str(i)+', left']['x'],y=method[1]['SC'+str(i)+', left']['y'])
                self.tele_r_ang = lambda i,t: pack.functions.get_tele_SS(False,False,i,t,'r',x=method[1]['SC'+str(i)+', right']['x'],y=method[1]['SC'+str(i)+', right']['y'])

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
            beam_l_ang=[]
            beam_r_ang=[]
            for i in range(1,4):
                t_l =  getattr(getattr(ret[0],'l'),'i'+str(i)).tele_ang
                t_r =  getattr(getattr(ret[0],'r'),'i'+str(i)).tele_ang
                b_r =  getattr(getattr(ret[0],'r'),'i'+str(i)).PAAM_ang
                b_r =  getattr(getattr(ret[0],'r'),'i'+str(i)).PAAM_ang
                tele_l_ang.append(methods.interpolate(t_l[0][0],t_l[0][1]))
                tele_r_ang.append(methods.interpolate(t_r[0][0],t_r[0][1]))
                beam_l_ang.append(methods.interpolate(t_l[0][0],t_l[0][1]))
                beam_r_ang.append(methods.interpolate(t_r[0][0],t_r[0][1]))

            self.tele_l_ang_func = lambda i,t: tele_l_ang[i-1](t)
            self.tele_r_ang_func = lambda i,t: tele_r_ang[i-1](t)
            self.beam_l_ang_func = lambda i,t: beam_l_ang[i-1](t)
            self.beam_r_ang_func = lambda i,t: beam_r_ang[i-1](t)

#        elif 'pointLISA.utils.Object' in str(type(self.aimset.inp)):
#            option = self.tele_control+'_'+self.PAAM_control+'__'+self.option_tele+'_'+self.option_PAAM
#            setting = self.tele_method_solve+'_'+self.PAAM_method_solve+'__'+self.optimize_PAAM+'_'+str(self.optimize_PAAM_value).replace('.','d')
#
#            self.inp = getattr(getattr(self.aimset.inp,option),setting)
#            print('Reading imported telescope angles')
#
#            XL = []
#            YL = []
#            XR = []
#            YR = []
#            
#            tele_l_ang_func = []
#            tele_r_ang_func = []
#            beam_l_ang_func = []
#            beam_r_ang_func = []
#
#            for i in range(1,4):
#                [xlt,ylt] = getattr(getattr(self.inp.l,'i'+str(i)),'tele_l_ang')[0]
#                [xrt,yrt] = getattr(getattr(self.inp.r,'i'+str(i)),'tele_r_ang')[0]
#                [xlb,ylb] = getattr(getattr(self.inp.l,'i'+str(i)),'beam_l_ang')[0]
#                [xrb,yrb] = getattr(getattr(self.inp.r,'i'+str(i)),'beam_r_ang')[0]
#
#                tele_l_ang_func.append(methods.interpolate(xlt,ylt))
#                tele_r_ang_func.append(methods.interpolate(xrt,yrt))
#                beam_l_ang_func.append(methods.interpolate(xlb,ylb))
#                beam_r_ang_func.append(methods.interpolate(xrb,yrb))
#
#            self.tele_l_ang = lambda i,t: tele_l_ang_func[i-1](t)
#            self.tele_r_ang = lambda i,t: tele_r_ang_func[i-1](t)
#            self.beam_l_ang = lambda i,t: beam_l_ang_func[i-1](t)
#            self.beam_r_ang = lambda i,t: beam_r_ang_func[i-1](t)            
#
        
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
                print("Sampling and fitting telescope angles")
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
        print('PAAM pointing strategy: '+option)

        if option=='center':
            #ang_l = lambda i,t: output.PAAM_center_calc(self.aim0,i,t,tele_l=self.tele_l_ang(i,t),tele_r=self.tele_r_ang(utils.i_slr(i)[1],t),lim=self.aimset.limits.yoff,method=self.aimset.PAAM_method_solve,para=self.aimset.optimize_PAAM,value=self.aimset.optimize_PAAM_value)[0][0]
            ang_l = lambda i,t: output.PAAM_center_calc(self.aim0,i,t,tele_l=self.tele_l_ang(i,t),tele_r=self.tele_r_ang(utils.i_slr(i)[1],t),lim=self.aimset.limit_yoff,method=self.aimset.PAAM_method_solve,para=self.aimset.optimize_PAAM,value=self.aimset.optimize_PAAM_value,margin=self.optimize_PAAM_margin)[0][0]
            #ang_r = lambda i,t: output.PAAM_center_calc(self.aim0,utils.i_slr(i)[2],t,tele_l=self.tele_l_ang(utils.i_slr(i)[2],t),tele_r=self.tele_r_ang(i,t),lim=self.aimset.limits.yoff,method=self.aimset.PAAM_method_solve,para=self.aimset.optimize_PAAM,value=self.aimset.optimize_PAAM_value)[0][1]
            ang_r = lambda i,t: output.PAAM_center_calc(self.aim0,utils.i_slr(i)[2],t,tele_l=self.tele_l_ang(utils.i_slr(i)[2],t),tele_r=self.tele_r_ang(i,t),lim=self.aimset.limit_yoff,method=self.aimset.PAAM_method_solve,para=self.aimset.optimize_PAAM,value=self.aimset.optimize_PAAM_value,margin=self.optimize_PAAM_margin)[0][1]
            
            #ang_l = lambda i,t: output.PAAM_center_calc(self.aim0,i,t,lim=self.aimset.limits.yoff)[0][0]
            #ang_r = lambda i,t: output.PAAM_center_calc(self.aim0,utils.i_slr(i)[2],t,lim=self.aimset.limits.yoff)[0][1]


        elif option=='wavefront':
            ang_l = lambda i,t: output.PAAM_wavefront_calc(self,i,t,'l',lim=self.aimset.limit_angy)
            ang_r = lambda i,t: output.PAAM_wavefront_calc(self,i,t,'r',lim=self.aimset.limit_angy)


            #ang_l = lambda i,t: output.PAAM_wavefront_calc(self,i,t,'l',lim=self.aimset.limit_angy)
            #ang_r = lambda i,t: output.PAAM_wavefront_calc(self,i,t,'r',lim=self.aimset.limit_angy)
            #ang_l = lambda i,t: methods.rotate_PAA_wavefront(self.data,self.aim_old,i,t,'l','angy')
            #ang_r = lambda i,t: methods.rotate_PAA_wavefront(self.data,self.aim_old,i,t,'r','angy')

        self.PAAM_option = option
        return [ang_l,ang_r]


    def PAAM_aim(self,method=False,dt=3600*24,jitter=False,tau=1,mode='overdamped',PAAM_ang_extra=False,option='center'):
        if method==False:
            method = self.PAAM_method
        else:
            self.PAAM_method = method

        print('The PAAM control method is: ' +method)
        print(' ')

        # Obtaining PAAM angles for 'fc' (full_control), 'nc' (no_control) and 'SS' (step and stair)
        
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
                #self.aim_sampled.get_coordinate_systems(iteration_val=False,option='self')
            else:
                self.aim_sampled=False

        #self.get_coordinate_systems(iteration_val=False,option='self')
        if self.inp!=False:
            self.aim_sampled = self
        return self

    def sample(self): 
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
        # Calculating new pointing vectors and coordinate system
        tele_l_coor = methods.coor_tele(self.data,i,t,tele_l_ang(i,t))
        tele_r_coor = methods.coor_tele(self.data,i,t,tele_r_ang(i,t))
        tele_l_vec = LA.unit(methods.coor_tele(self.data,i,t,tele_l_ang(i,t))[0])*self.data.L_tele
        tele_r_vec = LA.unit(methods.coor_tele(self.data,i,t,tele_r_ang(i,t))[0])*self.data.L_tele
        tele_l_start = tele_l_vec+np.array(self.data.putp(i,t))
        tele_r_start = tele_r_vec+np.array(self.data.putp(i,t))

        return [[tele_l_coor,tele_r_coor],[tele_l_vec,tele_r_vec],[tele_l_start,tele_r_start]]


    def get_beam_coor(self,i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang):
        # Calculating new pointing vectors and coordinate system
        beam_l_coor = methods.beam_coor_out(self.data,i,t,tele_l_ang(i,t),beam_l_ang(i,t),self.offset_tele['l'])
        beam_r_coor = methods.beam_coor_out(self.data,i,t,tele_r_ang(i,t),beam_r_ang(i,t),self.offset_tele['r'])

        # Calculating the Transmitted beam direction and position of the telescope aperture
        beam_l_direction = beam_l_coor[0]
        beam_r_direction = beam_r_coor[0]
        beam_l_start = beam_l_direction+np.array(self.data.putp(i,t))
        beam_r_start = beam_r_direction+np.array(self.data.putp(i,t))

        return [[beam_l_coor,beam_r_coor],[beam_l_direction,beam_r_direction],[beam_l_start,beam_r_start]]

    def get_coordinate_systems(self,iteration_val=False,option='self'):
        if iteration_val==False:
            tele_l_ang = self.tele_l_ang
            tele_r_ang = self.tele_r_ang
            beam_l_ang = self.beam_l_ang
            beam_r_ang = self.beam_r_ang
        else:
            [[tele_l_ang,tele_r_ang],[beam_l_ang,beam_r_ang]] = self.get_funcions_from_sampling(self.get_sampled_pointing(option=option))

        self.tele_l_coor = lambda i,t: self.get_tele_coor(i,t,tele_l_ang,tele_r_ang)[0][0]
        self.tele_r_coor = lambda i,t: self.get_tele_coor(i,t,tele_l_ang,tele_r_ang)[0][1]
        self.tele_l_vec = lambda i,t: self.get_tele_coor(i,t,tele_l_ang,tele_r_ang)[1][0]
        self.tele_r_vec = lambda i,t: self.get_tele_coor(i,t,tele_l_ang,tele_r_ang)[1][1]
        self.tele_l_start = lambda i,t: self.get_tele_coor(i,t,tele_l_ang,tele_r_ang)[2][0]
        self.tele_r_start = lambda i,t: self.get_tele_coor(i,t,tele_l_ang,tele_r_ang)[2][1]


        self.beam_l_coor = lambda i,t: self.get_beam_coor(i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang)[0][0]
        self.beam_r_coor = lambda i,t: self.get_beam_coor(i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang)[0][1]
        self.beam_l_direction = lambda i,t: self.get_beam_coor(i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang)[1][0]
        self.beam_r_direction = lambda i,t: self.get_beam_coor(i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang)[1][1]
        self.beam_l_start = lambda i,t: self.get_beam_coor(i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang)[2][0]
        self.beam_r_start = lambda i,t: self.get_beam_coor(i,t,tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang)[2][1]

        self.tele_l_ang_calc = tele_l_ang
        self.tele_r_ang_calc = tele_r_ang
        self.beam_l_ang_calc = beam_l_ang
        self.beam_r_ang_calc = beam_r_ang

        return 0 #[tele_l_ang,tele_r_ang,beam_l_ang,beam_r_ang,tele_l_coor,tele_r_coor,tele_l_vec,tele_r_vec,beam_l_coor,beam_r_coor,beam_l_direction,beam_r_direction,beam_l_start,beam_r_start]


    def get_aim_accuracy(self,i,t,side,component=False,option='wavefront'):
        [i_self,i_left,i_right] = utils.i_slr(i)
        if option=='center':
            if component=='tele':
                if side=='l':
                    tdel=self.data.L_rl_func_tot(i_self,t)
                    if self.data.calc_method=='Waluschka':
                        tdel0=tdel
                    elif self.data.calc_method=='Abram':
                        tdel0=0
                    end = self.aim_old.tele_l_start(i_self,t-tdel0)
                    start = self.aim_old.tele_r_start(i_left,t-tdel)
                    coor_start = self.aim_old.beam_r_coor(i_left,t-tdel)
                    coor_end = self.aim_old.tele_l_coor(i_self,t)
                elif side=='r':
                    tdel = self.data.L_rr_func_tot(i_self,t)
                    if self.data.calc_method=='Waluschka':
                        tdel0=tdel
                    elif self.data.calc_method=='Abram':
                        tdel0=0
                    end = self.aim_old.tele_r_start(i_self,t-tdel0)
                    start = self.aim_old.tele_l_start(i_right,t-tdel)
                    coor_start = self.aim_old.beam_l_coor(i_right,t-tdel)
                    coor_end = self.aim_old.tele_r_coor(i_self,t)

                [z,y,x] = LA.matmul(coor_start,end-start)
                angx = np.sign(x)*abs(np.arctan(x/z))
                angy = np.sign(y)*abs(np.arctan(y/z))
                delay = tdel
            elif component=='PAAM':
                if side=='l':
                    tdel = self.data.L_sr_func_tot(i_left,t)
                    if self.data.calc_method=='Waluschka':
                        tdel0=tdel
                    elif self.data.calc_method=='Abram':
                        tdel0=0
                    start = self.aim_old.tele_r_start(i_left,t+tdel0)
                    end = self.aim_old.tele_l_start(i_self,t+tdel)
                    coor_start = self.aim_old.beam_r_coor(i_left,t)
                elif side=='r':
                    tdel = self.data.L_sl_func_tot(i_right,t)
                    if self.data.calc_method=='Waluschka':
                        tdel0=tdel
                    elif self.data.calc_method=='Abram':
                        tdel0=0
                    start = self.aim_old.tele_l_start(i_right,t+tdel0)
                    end = self.aim_old.tele_r_start(i_self,t+tdel)
                    coor_start = self.aim_old.beam_l_coor(i_right,t)

                [z,y,x] = LA.matmul(coor_start,end-start)
                angx = np.sign(x)*abs(np.arctan(x/z))
                angy = np.sign(y)*abs(np.arctan(y/z))
                delay = tdel

            return [angx,-angy,delay]

        elif option=='wavefront':
            if component=='tele':
                if side=='l':
                    tdel = self.data.L_rl_func_tot(i_self,t)
                    if self.data.calc_method=='Waluschka':
                        tdel0=tdel
                    elif self.data.calc_method=='Abram':
                        tdel0=0
                    coor_start = self.aim_old.tele_r_coor(i_left,t-tdel)
                    coor_end = self.aim_old.tele_l_coor(i_self,t)
                    direct = self.aim_old.beam_r_direction(i_left,t-tdel)
                    start = self.aim_old.tele_r_start(i_left,t-tdel)
                    end = self.aim_old.tele_l_start(i_self,t-tdel0)
                elif side=='r':
                    tdel = self.data.L_rr_func_tot(i_self,t)
                    if self.data.calc_method=='Waluschka':
                        tdel0=tdel
                    elif self.data.calc_method=='Abram':
                        tdel0=0
                    coor_start = self.aim_old.tele_l_coor(i_right,t-tdel)
                    coor_end = self.aim_old.tele_r_coor(i_self,t)
                    direct = self.aim_old.beam_l_direction(i_right,t-tdel)
                    start = self.aim_old.tele_l_start(i_right,t-tdel)
                    end = self.aim_old.tele_r_start(i_self,t-tdel0)

                [zoff,yoff,xoff] = LA.matmul(coor_start,end-start)
                R = zoff #Not precise
                R_vec = np.array([(R**2-xoff**2-yoff**2)**0.5,yoff,xoff])
                R_vec_origin = LA.matmul(np.linalg.inv(coor_start),R_vec)
                [z,y,x] = LA.matmul(coor_end,-R_vec_origin)
                #[z,y,x]=LA.matmul(coor_end,-direct)
                angx = np.sign(x)*abs(np.arctan(x/z))
                angy = np.sign(y)*abs(np.arctan(y/z))
                delay=tdel

            elif component=='PAAM':
                angy_solve = lambda PAAM_ang: methods.get_wavefront_parallel(self.wfe,self.aim_old,i_self,t,side,PAAM_ang,'angy',mode='opposite')
                try:
                    angy =  scipy.optimize.brentq(angy_solve,-1e-1,1e-1)
                except ValueError:
                    angy=np.nan
                delay=False
                angx=False

        return [angx,angy,delay]
    
#    def copy_aim(self,aim_old,option=False): #...copying all essential parameters
#        for attr in aim_old.__dict__.keys():
#
#            if attr not in self.__dict__.keys() or attr=='wfe':
#                t = str(aim_old.__dict__[attr])
#                go=False
#                if option!=False:
#                    if option=='new_angles':
#                        if 'instance' in t:
#                            print(t)
#                            if 'WFE' in t:
#                                if type(self.wfe)==bool:
#                                    go=True
#                                else:
#                                    go=False
#                            else:
#                                go=False
#                        elif 'function' in t:
#                            go =False
#                        else:
#                            go=True
#                else:
#                    go=True
#                if go==True:
#                    #print(attr)
#                    setattr(self,attr,aim_old.__dict__[attr])
#
#    return 0

