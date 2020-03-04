from imports import * 
import inspect
import numpy as np
import methods
from pointLISA import * 
import scipy
import utils
class OUTPUT():
    '''This class uses an AIM object and calculate some useful properties (output) as attributes of an OUTPUT object'''
    def __init__(self,aim=False,**kwargs):
        from pointLISA import * 
        import pointLISA.utils as utils
        import pointLISA

        globals()['utils'] = pointLISA.utils
        if aim!=False:
            self.aim = aim
        else:
            self.aim = utils.Object()
            self.aim.data = kwargs['data']

        ### Obtain parameters
        for k, value in parameters.__dict__.items():
            if '__' not in k:
                globals()[k] = value
                setattr(self,k,value)
    
    def pupil(self,Nbins=2): #Aperture
        '''Creates pixels on the aperture which can be useful to examine the variation of a property over the aperture'''
        try:
            del_at = ['xlist','ylist','Deltax','Deltay','Nbinsx','Nbinsy']
            for at in del_at: 
                delattr(self,at)
        except:
            pass

        D_calc=self.aim.data.D

        xlist = np.linspace(-D_calc*0.5,D_calc*0.5,Nbins+1)
        self.xlist = xlist[0:-1]+0.5*(xlist[1]-xlist[0])
        self.ylist = self.xlist
        self.Deltax = self.xlist[1]-self.xlist[0]
        self.Deltay = self.ylist[1]-self.ylist[0]
        self.Nbinsx = len(self.xlist)
        self.Nbinsy = len(self.ylist)

    def w(self,z):
        '''Beamwaist as a function of z (z=coordinate along beamline)'''
        zR = np.pi*(self.w0_laser**2)/self.labda

        return self.w0_laser*((1+((z/zR)**2))**0.5)

    def R(self,z,guess=False):
        '''The radius of curvasture R as a function of z'''
        if z!=np.nan:
            zR = np.pi*(self.w0_laser**2)/self.labda

            if guess==False:
                return abs(z*(1+((zR/z)**2)))

            elif guess==True:
                return z
        else:
            return np.nan
       
    def u0(self,ksi):#...for gaussian beam, adjust when implementing other beams
        '''Beam profile (amplitude)'''
        w = self.w(0.0)
        [x,y] = ksi

        return np.exp(-((x**2+y**2)/(w**2)))

    def aperture(self,xlist,ylist,function,dType=np.float64): # Creates matrix of function over an aperture (circle)
        '''Obtains the values of a parameter (function) over the entire aperture'''
        if type(xlist)==bool:
            if xlist==False:
                xlist = self.xlist
        if type(ylist)==bool:
            if ylist==False:
                ylist = self.ylist

        Nbins = len(xlist)
        ps = np.empty((Nbins,Nbins),dtype=dType)
        for i in range(0,len(xlist)):
            for j in range(0,len(ylist)):
                x = xlist[i]
                y = ylist[j]
                if x**2+y**2<=0.25*(self.D**2):
                    val = function(x,y)
                    if type(val)==list:
                        val=val[0]
                    ps[i,j] = val

                else:
                    ps[i,j] = np.nan

        return ps
    
    def z_solve(self,x,y,z,calc_R=False,ret='piston',R_guess=True):
        '''Solves the photon traveling distance from the point of transmitting to the point [x,y,z]'''
        try:
            if z!=np.nan:
                x = np.float64(x)
                y = np.float64(y)
                z = np.float64(z)
                f_solve = lambda dz: (self.R(z+dz,guess=R_guess) - (self.R(z+dz,guess=R_guess)**2 - (x**2+y**2))**0.5) - dz
                f_solve_2 = lambda dz: (z- (((z+dz)**2 - x**2 -y**2 )**0.5))
                dz_sol = scipy.optimize.brentq(f_solve,-0.5*z,0.5*z,xtol=1e-64)
            else:
                dz_sol=np.nan
                dz_sol_2=np.nan
                dz_sol_3=np.nan
                raise ValueError

            if calc_R==True:
                return self.R(z+dz_sol,guess=R_guess)
            else:
                if ret=='piston':
                    return z+dz_sol
                elif ret=='all':
                    return [z+dz_sol,dz_sol]

        except RuntimeError:
            if ret=='piston':
                return np.nan
            elif ret=='all':
                return [np.nan,np.nan]


    def get_coordinate_systems(self,speed=0):
        '''Gets the coordinate systems'''
        if speed==0:
            self.aim.get_coordinate_systems()

        elif speed>0:
           if self.aim.sampled==False:
               self.aim.aim_sampled = self.aim.sample()
        try:
            self.aim.aim_sampled.get_coordinate_systems()
        except:
            self.aim.get_coordinate_systems()
        return self.aim

    def add_attribute(self,e,pos):
        ''''''
        func = 'get_'+str(e).split("attribute '")[-1].replace("'","")
        pos = getattr(self,func)(pos)

        return pos

    ### All 'get_...' functions defined below will obtain an output value
     
    def get_end_ksi(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.ksi[0]*pos.coor_end[2]+pos.ksi[1]*pos.coor_end[1]
                print('ret:')
                print(ret)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_offset(self,pos): #Done
        check=False
        while check==False:
            try:
                if pos.side=='l':
                    ret = pos.offset_l
                elif pos.side=='r':
                    ret = pos.offset_r
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_off(self,pos): #Done
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_startbeam__send,pos.arm__send)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                setattr(pos,'xoff',ret[2])
                setattr(pos,'yoff',ret[1])
                setattr(pos,'zoff',ret[0])
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos
   
    def get_tele_ang(self,pos): #Done
        check=False
        while check==False:
            try:
                if pos.side=='l':
                    ret = self.aim.tele_l_ang(pos.i_self,pos.t)
                elif pos.side=='r':
                    ret = self.aim.tele_r_ang(pos.i_self,pos.t)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_PAAM_ang(self,pos): #Done
        check=False
        while check==False:
            try:
                if pos.side=='l':
                    ret = self.aim.beam_l_ang(pos.i_self,pos.t)
                elif pos.side=='r':
                    ret = self.aim.beam_r_ang(pos.i_self,pos.t)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos
    
    def get_invside(self,pos): #Done
        check=False
        while check==False:
            try:
                if pos.side=='l':
                    ret='r'
                elif pos.side=='r':
                    ret='l'
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_i_opp(self,pos): #Done
        check=False
        while check==False:
            try:
                if pos.side=='l':
                    ret=pos.i_left
                elif pos.side=='r':
                    ret=pos.i_right
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_coor_starttele(self,pos):
        check=False
        while check==False:
            try:
                if pos.mode=='send':
                    ret=get_coor_tele(pos.aim,pos.i_self,pos.t,pos.side,tele_angle=pos.tele_angle_start)
                elif pos.mode=='rec':
                    ret=get_coor_tele(pos.aim,pos.i_opp,pos.t-pos.tdel,pos.invside,tele_angle=pos.tele_angle_start)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_coor_startbeam__send(self,pos):
        check=False
        while check==False:
            try:
                if pos.mode=='send':
                    ret=get_coor_beam_out__send(pos.aim,pos.i_send,pos.t,pos.side,tele_angle=pos.tele_angle_start,beam_angle=pos.beam_angle_start,offset=pos.offset_start)
                elif pos.mode=='rec':
                    ret=get_coor_beam_out__send(pos.aim,pos.i_send,pos.t-pos.tdel,pos.invside,tele_angle=pos.tele_angle_start,beam_angle=pos.beam_angle_start,offset=pos.offset_start)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_vec_startbeam__send(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.coor_startbeam__send[0]
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_vec_startbeam__sun(self,pos):
        check=False
        while check==False:
            try:
                v = pos.coor_startbeam__send
                if pos.mode=='send':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t+pos.tdel0,v,reverse=False)
                elif pos.mode=='rec':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t-pos.tdel,v,reverse=False)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_vec_startbeam__rec(self,pos):
        check=False
        while check==False:
            try:
                v = pos.coor_startbeam__sun
                if pos.mode=='send':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t+pos.tdel,v,reverse=True)
                elif pos.mode=='rec':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t-pos.tdel0,v,reverse=True)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos
    
    def get_coor_endbeam__send(self,pos):
        check=False
        while check==False:
            try:
                if pos.mode=='send':
                    ret=get_coor_beam_out__send(pos.aim,pos.i_send,pos.t+pos.tdel,pos.side,tele_angle=pos.tele_angle_start,beam_angle=pos.beam_angle_start,offset=pos.offset_start)
                elif pos.mode=='rec':
                    ret = get_coor_beam_out__send(pos.aim,pos.i_send,pos.t-pos.tdel0,pos.invside,tele_angle=pos.tele_angle_start,beam_angle=pos.beam_angle_start,offset=pos.offset_start)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_vec_endbeam__send(self,pos): #send frame, Sun CS
        check=False
        while check==False:
            try:
                ret = pos.coor_endbeam__send[0]
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_vec_endbeam__sun(self,pos): #Sun frame, Sun CS
        check=False
        while check==False:
            try:
                v = pos.coor_endbeam__send
                
                if pos.mode=='send':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t+pos.tdel0,v,reverse=False)
                elif pos.mode=='rec':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t-pos.tdel,v,reverse=False)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_vec_endbeam__rec(self,pos): #Rec frame, Sun CS
        check=False
        while check==False:
            try:
                v = pos.coor_startbeam__sun
                if pos.mode=='send':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t+pos.tdel,v,reverse=True)
                elif pos.mode=='rec':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t-pos.tdel0,v,reverse=True)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_coor_end(self,pos):
        check=False
        while check==False:
            try:
                if pos.mode=='send':
                    ret = get_coor_tele(pos.aim,pos.i_opp,pos.t+pos.tdel,pos.invside,tele_angle=pos.tele_angle_end)
                elif pos.mode=='rec':
                    ret = get_coor_tele(pos.aim,pos.i_self,pos.t-pos.tdel0,pos.side,tele_angle=pos.tele_angle_end)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_start(self,pos):
        check=False
        while check==False:
            try:
                if pos.mode=='send':
                    ret = get_start_calc(pos.aim,pos.i_self,pos.t+pos.tdel0,pos.side,pos.tele_angle_start)
                elif pos.mode=='rec':
                    ret = get_start_calc(pos.aim,pos.i_opp,pos.t-pos.tdel,pos.invside,pos.tele_angle_start)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_end(self,pos):
        check=False
        while check==False:
            try:
                if pos.mode=='send':
                    ret = get_start_calc(pos.aim,pos.i_opp,pos.t+pos.tdel,pos.invside,pos.tele_angle_end)
                elif pos.mode=='rec':
                    ret = get_start_calc(pos.aim,pos.i_self,pos.t-pos.tdel0,pos.side,pos.tele_angle_end)
                if pos.ksi!=[0,0]:
                    ret = ret+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_xoff(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.off[2]
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_yoff(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.off[1]
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_zoff(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.off[0]
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_startend__sun(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.end-pos.start
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                self.add_attribute(e,pos)
        return pos

    def get_startend__send(self,pos):
        check=False
        while check==False:
            try:
                v = pos.startend__sun
                if pos.mode=='send':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t+pos.tdel0,v,reverse=True)
                elif pos.mode=='rec':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t-pos.tdel,v,reverse=True)            
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_startend__rec(self,pos):
        check=False
        while check==False:
            try:
                v = pos.startend__sun
                if pos.mode=='send':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t+pos.tdel,v,reverse=True)
                elif pos.mode=='rec':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t-pos.tdel0,v,reverse=True)            
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

#...hier gebleven










    def get_r(self,pos):
        check=False
        while check==False:
            try:
                ret = (pos.xoff**2+pos.yoff**2)**0.5
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos


    def get_arm__rec(self,pos):
        check=False
        while check==False:
            try:
                if pos.mode=='send':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t+pos.tdel,pos.startend__sun,reverse=True)
                elif pos.mode=='rec':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t-pos.tdel0,pos.startend__sun,reverse=True)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos
   
    def get_arm__send(self,pos):
        check=False
        while check==False:
            try:
                if pos.mode=='send':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t+pos.tdel0,pos.startend__sun,reverse=True) #.......check on time it emits!!!
                elif pos.mode=='rec':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t-pos.tdel,pos.startend__sun,reverse=True)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                self.add_attribute(e,pos)
        return pos


    def get_arm_tele_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_end,-pos.arm__rec)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos
    
    def get_arm_tele_send(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_starttele,pos.arm__send)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angx_arm_tele_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = np.arctan(pos.arm_tele_rec[2]/pos.arm_tele_rec[0])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angy_arm_tele_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = np.arctan(pos.arm_tele_rec[1]/pos.arm_tele_rec[0])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos
    
    def get_ang_arm_tele_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = np.arctan(((pos.arm_tele_rec[2]**2+pos.arm_tele_rec[1]**2)**0.5)/pos.arm_tele_rec[0])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angx_arm_tele_send(self,pos):
        check=False
        while check==False:
            try:
                ret = np.arctan(pos.arm_tele_send[2]/pos.arm_tele_send[0])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angy_arm_tele_send(self,pos):
        check=False
        while check==False:
            try:
                ret = np.arctan(pos.arm_tele_send[1]/pos.arm_tele_send[0])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos















    def get_R(self,pos,precision=0):
        check=False
        while check==False:
            try:
                if precision==0:
                    ret=pos.zoff
                elif precision==1:
                    ret = self.R(pos.piston)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos) 
        return pos

    def get_R_vec_beam__send(self,pos): # in Sun coordinate system and send inertial frame!!!
        check=False
        while check==False:
            try:
                vec = np.array([(pos.R**2-pos.xoff**2-pos.yoff**2)**0.5,pos.yoff,pos.xoff]) # In beam frame
                ret = LA.matmul(np.linalg.inv(pos.coor_startbeam__send),vec) 
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_R_vec_beam__sun(self,pos):
        check=False
        while check==False:
            try:
                v = pos.R_vec_beam__send
                if pos.mode=='send':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t+pos.tdel0,v,reverse=False)
                elif pos.mode=='rec':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t-pos.tdel,v,reverse=False)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_R_vec_beam__rec(self,pos):
        check=False
        while check==False:
            try:
                v = pos.R_vec_beam__sun
                if pos.mode=='send':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t+pos.tdel,v,reverse=True)
                elif pos.mode=='rec':
                    ret = methods.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t-pos.tdel0,v,reverse=True)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos


    def get_R_vec_tele_send(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_starttele,R_vec_beam__send)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos
     

#    def get_R_vec_origin(self,pos): #...look
#        check=False
#        while check==False:
#            try:
#                ret = LA.matmul(np.linalg.inv(pos.coor_startbeam_out),pos.R_vec_beam_send)
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos

    def get_R_vec_tele_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_end,pos.R_vec_beam__rec)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angx_R_vec_tele_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = np.arctan(pos.R_vec_tele_rec[2]/pos.R_vec_tele_rec[0])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angy_R_vec_tele_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = np.arctan(pos.R_vec_tele_rec[1]/pos.R_vec_tele_rec[0])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_ang_R_vec_tele_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = np.arctan(((pos.R_vec_tele_rec[2]**2+pos.R_vec_tele_rec[1]**2)**0.5)/pos.R_vec_tele_rec[0])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos









    def get_beam_receive_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_end,-pos.vec_endbeam__rec)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angx_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.beam_receive_rec[2])*abs(np.arctan(pos.beam_receive_rec[2]/pos.beam_receive_rec[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angy_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.beam_receive_rec[1])*abs(np.arctan(pos.beam_receive_rec[1]/pos.beam_receive_rec[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_beam_receive_send(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_starttele,pos.vec_endbeam__send)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos
    
    def get_angx_send(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.beam_receive_send[2])*abs(np.arctan(pos.beam_receive_send[2]/pos.beam_receive_send[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angy_send(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.beam_receive_send[1])*abs(np.arctan(pos.beam_receive_send[1]/pos.beam_receive_send[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

#    def get_R_vec_tele_rec(self,pos): #This vector is reversed (pointed away from receiving telescope)
#        check=False
#        while check==False:
#            try:
#                ret = LA.matmul(pos.coor_end,-pos.R_vec_origin)
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos
#
#    def get_R_vec_tele_rec_ab(self,pos): #This vector is reversed (pointed away from receiving telescope)
#        check=False
#        while check==False:
#            try:
#                ret = aberration(pos,pos.R_vec_tele_rec,mode='rec')
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos

  
    def get_tele_vec(self,pos): #This vector is reversed (pointed away from receiving telescope)
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_starttele,-pos.coor_end[0])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos
#
#    def get_tele_vec_ab(self,pos): #This vector is reversed (pointed away from receiving telescope)
#        check=False
#        while check==False:
#            try:
#                ret = aberration(pos,pos.tele_vec,mode='send')
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos

    def get_R_vec_beam_send__send(self,pos): # in send coordinate system and send inertial frame!!!
        check=False
        while check==False:
            try:
                v_sun = pos.R_vec_beam__send
                if pos.mode=='send':
                    v_send = methods.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t+pos.tdel0,v_sun,reverse=True)
                elif pos.mode=='rec':
                    v_send = methods.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t-pos.tdel,v_sun,reverse=True)
                ret = LA.matmul(pos.coor_startbeam__send,v_send)        
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angx_R(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.R_vec_beam_send__send[2])*abs(np.arctan(pos.R_vec_beam_send__send[2]/pos.R_vec_beam_send__send[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angy_R(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.R_vec_beam_send__send[1])*abs(np.arctan(pos.R_vec_beam_send__send[1]/pos.R_vec_beam_send__send[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos
    
#    def get_R_r2_wf_rec(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = pos.R_vec_tele_rec[1]**2+pos.R_vec_tele_rec[2]**2
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos

    def get_angx_wf_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.R_vec_tele_rec[2])*abs(np.arctan(pos.R_vec_tele_rec[2]/pos.R_vec_tele_rec[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angy_wf_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.R_vec_tele_rec[1])*abs(np.arctan(pos.R_vec_tele_rec[1]/pos.R_vec_tele_rec[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angx_wf_send(self,pos):
        check=False
        while check==False:
            try:
                endtele = -pos.coor_endtele[0]
                endtele_send = LA.matmul(pos.coor_starttele,endtele)
                R = pos.R_vec_beam_send__send
                angx_tele = np.sign(endtele_send[2])*np.arctan(abs(endtele_send[2]/endtele_send[0]))
                angx_R = np.sign(R[2])*np.arctan(abs(R[2]/R[0]))
                
                ret = angx_R - angx_tele
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angy_wf_send(self,pos):
        check=False
        while check==False:
            try:
                endtele = -pos.coor_endtele[0]
                endtele_send = LA.matmul(pos.coor_starttele,endtele)
                R = pos.R_vec_beam_send__send
                angx_tele = np.sign(endtele_send[1])*np.arctan(abs(endtele_send[1]/endtele_send[0]))
                angx_R = np.sign(R[1])*np.arctan(abs(R[1]/R[0]))

                ret = angx_R - angx_tele
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos
  
#    def get_angx_wf_send(self,pos): #...check
#        check=False
#        while check==False:
#            try:
#                tele_end = LA.matmul(pos.coor_starttele,-pos.coor_end[0]) # tele_rec is reversed
#                tele_end_ab = aberration(pos,tele_end,mode='send')
#                angx_tele_rec = np.sign(tele_end_ab[2])*np.arctan(abs(tele_end_ab[2]/tele_end_ab[0]))
#                angx_R_vec_tele_send = np.sign(pos.R_vec_tele_send_ab[2])*np.arctan(abs(pos.R_vec_tele_send_ab[2]/pos.R_vec_tele_send_ab[0]))
#                ret = angx_tele_rec - angx_R_vec_tele_send
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos
#
#    def get_angy_wf_send(self,pos):
#        check=False
#        while check==False:
#            try:
#                tele_end = LA.matmul(pos.coor_starttele,-pos.coor_end[0]) # tele_rec is reversed
#                tele_end_ab = aberration(pos,tele_end,mode='send')
#                angy_tele_rec = np.sign(tele_end_ab[1])*np.arctan(abs(tele_end_ab[1]/tele_end_ab[0]))
#                angy_R_vec_tele_send = np.sign(pos.R_vec_tele_send_ab[1])*np.arctan(abs(pos.R_vec_tele_send_ab[1]/pos.R_vec_tele_send_ab[0]))
#                ret = angy_tele_rec - angy_R_vec_tele_send
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos
#
#    def get_angx_wf_rec(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = np.sign(pos.R_vec_tele_rec_ab[2])*np.arctan(abs(pos.R_vec_tele_rec_ab[2]/pos.R_vec_tele_rec_ab[0]))
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos
#
#    def get_angy_wf_rec(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = np.sign(pos.R_vec_tele_rec_ab[1])*np.arctan(abs(pos.R_vec_tele_rec_ab[1]/pos.R_vec_tele_rec_ab[0]))
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #bprint(e)
#                self.add_attribute(e,pos)
#        return pos

### Hier gebleven
#    def get_beam_direction_origin(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = np.array(pos.coor_startbeam[0])
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos
#
#    def get_beam_direction_rec(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = LA.matmul(pos.coor_end,pos.beam_direction_origin)
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos
    
    def get_pistonandz(self,pos):
        check=False
        while check==False:
            try:
                ret = self.z_solve(pos.xoff,pos.yoff,pos.zoff,ret='all')
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                self.add_attribute(e,pos)
        return pos

    def get_piston(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.pistonandz[0]
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos 

    def get_z_extra(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.pistonandz[1]
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

#    def get_FOV_beamline(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = np.arccos(-pos.beam_direction_rec[0]/np.linalg.norm(pos.beam_direction_rec))
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos

    def get_alpha(self,pos):
        check=False
        while check==False:
            try:
                ret = (abs(pos.angx_wf_rec)**2 +abs(pos.angy_wf_rec)**2)**0.5
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

#    def get_FOV_position(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = LA.angle(pos.start-pos.end,pos.coor_end[0])
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos
#
#    def get_tilt(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = pos.FOV_wavefront
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos
   
    def get_waist(self,pos):
        check=False
        while check==False:
            try:
                ret = self.w(pos.zoff)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos
    
    def get_FOVlim(self,pos):
        check=False
        while check==False:
            try:
                if pos.alpha>pos.aim.data.FOV:#pos.waist>pos.aim.data.FOV:
                    ret=0
                else:
                    ret=1
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_Ival(self,pos):
        check=False
        while check==False:
            try:
                ret = (pos.I0*np.exp((-2*(pos.xoff**2+pos.yoff**2))/(pos.waist**2)))*(pos.aim.data.w0_laser/pos.waist)

                #ret = (abs(pos.u)**2)[0]#*np.cos(pos.angx_rec)*np.cos(pos.angy_rec)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_Ivalx(self,pos):
        check=False
        while check==False:
            try:
                ret = (pos.I0*np.exp((-2*(pos.xoff**2))/(pos.waist**2)))*(pos.aim.data.w0_laser/pos.waist)

                #ret = (abs(pos.u)**2)[0]#*np.cos(pos.angx_rec)*np.cos(pos.angy_rec)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos


    
    def get_I(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.Ival*pos.FOVlim
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_I0(self,pos):
        check=False
        while check==False:
            try:
                ret = (self.P_L*np.pi*(self.w0_laser**2))/2.0
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_power(self,pos):
        check=False
        while check==False:
            try:
                try:
                    xlist=self.xlist
                    ylist=self.ylist
                except AttributeError:
                    self.pupil(Nbins=self.Nbins)
                    xlist=self.xlist
                    ylist=self.ylist

                if len(xlist)==1 and len(ylist)==1:
                    dksi = (self.D**2)*(np.pi/4.0)
                else:
                    dksi = (xlist[1]-xlist[0])*(ylist[1]-ylist[0])
                ret = dksi*pos.I
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_tdel(self,pos):
        check=False
        while check==False:
            try:
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
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_tdel0(self,pos):
        check=False
        while check==False:
            try:
                if pos.calc_method=='Abram':
                    ret=0
                elif pos.calc_method=='Waluschka':
                    ret = pos.tdel
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

#    def get_coor_starttele(self,pos):
#        check=False
#        while check==False:
#            try:
#                if pos.mode=='send':
#                    if pos.side=='l':
#                        ret = get_coor_tele(pos.aim,pos.i_self,pos.t+pos.tdel0,'l',tele_angle=pos.tele_angle_l)
#                    elif pos.side=='r':
#                        ret = get_coor_tele(pos.aim,pos.i_self,pos.t+pos.tdel0,'r',tele_angle=pos.tele_angle_r)
#                elif pos.mode=='rec':
#                    if pos.side=='l':
#                        ret = get_coor_tele(pos.aim,pos.i_left,pos.t-pos.tdel,'r',tele_angle=pos.tele_angle_r)
#                    elif pos.side=='r':
#                        ret = get_coor_tele(pos.aim,pos.i_right,pos.t-pos.tdel,'l',tele_angle=pos.tele_angle_l)
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError, e:
#                self.add_attribute(e,pos)
#        return pos
#
#    def get_coor_startbeam_out(self,pos): #...weghalen?
#        check=False
#        while check==False:
#            try:
#                if pos.mode=='send':
#                    if pos.side=='l':
#                        ret = get_coor_beam_out(pos.aim,pos.i_self,pos.t+pos.tdel0,'l',tele_angle=pos.tele_angle_l,beam_angle=pos.beam_angle_l,offset=pos.offset_l)
#                        
#                    elif pos.side=='r':
#                        ret = get_coor_beam_out(pos.aim,pos.i_self,pos.t+pos.tdel0,'r',tele_angle=pos.tele_angle_r,beam_angle=pos.beam_angle_r,offset=pos.offset_r)
#                elif pos.mode=='rec':
#                    if pos.side=='l':
#                        ret = get_coor_beam_out(pos.aim,pos.i_left,pos.t-pos.tdel,'r',tele_angle=pos.tele_angle_r,beam_angle=pos.beam_angle_r,offset=pos.offset_r)
#                    elif pos.side=='r':
#                        ret = get_coor_beam_out(pos.aim,pos.i_right,pos.t-pos.tdel,'l',tele_angle=pos.tele_angle_l,beam_angle=pos.beam_angle_l,offset=pos.offset_l)
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError, e:
#                self.add_attribute(e,pos)
#        return pos
#
    def get_coor_endtele(self,pos):
        check=False
        while check==False:
            try:
                if pos.mode=='send':
                    if pos.side=='l':
                        ret = get_coor_tele(pos.aim,pos.i_left,pos.t+pos.tdel,'r',tele_angle=pos.tele_angle_r)
                    elif pos.side=='r':
                        ret = get_coor_tele(pos.aim,pos.i_right,pos.t+pos.tdel,'l',tele_angle=pos.tele_angle_l)
                elif pos.mode=='rec':
                    if pos.side=='l':
                        ret = get_coor_tele(pos.aim,pos.i_self,pos.t-pos.tdel0,'l',tele_angle=pos.tele_angle_l)
                    elif pos.side=='r':
                        ret = get_coor_tele(pos.aim,pos.i_self,pos.t-pos.tdel0,'r',tele_angle=pos.tele_angle_r)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos
#
#    def get_start(self,pos):
#        check=False
#        while check==False:
#            try:
#                if pos.mode=='send':
#                    if pos.side=='l':
#                        ret = get_start_calc(pos.aim,pos.i_self,pos.t+pos.tdel0,'l',pos.tele_angle_l)
#                    elif pos.side=='r':
#                        ret = get_start_calc(pos.aim,pos.i_self,pos.t+pos.tdel0,'r',pos.tele_angle_r)
#                elif pos.mode=='rec':
#                    if pos.side=='l':
#                        ret = get_start_calc(pos.aim,pos.i_left,pos.t-pos.tdel,'r',pos.tele_angle_r)
#                    elif pos.side=='r':
#                        ret = get_start_calc(pos.aim,pos.i_right,pos.t-pos.tdel,'l',pos.tele_angle_l)
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError, e:
#                self.add_attribute(e,pos)
#        return pos
#
#    def get_end(self,pos):
#        check=False
#        while check==False:
#            try:
#                if pos.mode=='send':
#                    if pos.side=='l':
#                        ret = get_start_calc(pos.aim,pos.i_left,pos.t+pos.tdel,'r',pos.tele_angle_r)+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
#                    elif pos.side=='r':
#                        ret = get_start_calc(pos.aim,pos.i_right,pos.t+pos.tdel,'l',pos.tele_angle_l)+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
#                elif pos.mode=='rec':
#                    if pos.side=='l':
#                        ret = get_start_calc(pos.aim,pos.i_self,pos.t-pos.tdel0,'l',pos.tele_angle_l)+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
#                    elif pos.side=='r':
#                        ret = get_start_calc(pos.aim,pos.i_self,pos.t-pos.tdel0,'r',pos.tele_angle_r)+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError, e:
#                self.add_attribute(e,pos)
#        return pos
#
#    def get_direction(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = pos.coor_startbeam[0]
#                check=True
#            except AttributeError, e:
#                self.add_attribute(e,pos)
#        return pos


    def mean_var(self,i,t,side,ret,mode='mean',Nbins=False,tele_angle_l=False,tele_angle_r=False,beam_angle_l=False,beam_angle_r=False):
        '''Returns the output value for center (value at the center of the aperture), mean (mean value over the aperture), var (variance of the value over the aperture, mean_var (the mean and variance over the aperture) or mean_surface (matrix of values which represents the aperture pixels)'''
        if Nbins!=False:
            self.pupil(Nbins=Nbins)
        else:
            try:
                self.xlist
            except AttributeError:
                self.pupil()
        if type(ret)!=list:
            ret=[ret]

        func = lambda x,y: values(self,i,t,side,ret=ret,ksi=[x,y],tele_angle_l=tele_angle_l,tele_angle_r=tele_angle_r,beam_angle_l=beam_angle_l,beam_angle_r=beam_angle_r)

        if mode=='center':
            return getattr(func(0,0),ret[0])
        elif 'var' in mode or 'mean' in mode:
            A=[]
            for x in self.xlist:
                for y in self.ylist:
                    if ((x**2+y**2)**0.5)>=self.aim.data.D/2.0:
                        A.append(getattr(func(x,y),ret[0])*np.nan)
                    else:
                        A.append(getattr(func(x,y),ret[0]))
            if mode=='mean':
                return np.nanmean(A)
            elif mode=='var':
                return np.nanvar(A)/(np.float64(len(A) - A.count(np.nan)))
            elif mode=='mean_var':
                return np.array([np.nanmean(A),np.nanvar(A)/(np.float64(len(A) - A.count(np.nan)))])
            elif mode=='mean_surface':
                return A

### Write and save functions/values
    
    def t_calc(self,calc=False,**kwargs):
        '''Time array'''
        if calc==True:
            t0= kwargs.pop('t0',False)
            tend= kwargs.pop('tend',False)
            dt= kwargs.pop('dt',False)
            n= kwargs.pop('n',False)

            if dt==False:
                try:
                    dt = self.dt
                except:
                    dt = self.aim.data.t_all[1]-self.aim.data.t_all[0]
            if n!=False:
                tend = dt*n

            if t0==False:
                t0 = self.aim.data.t_all[3]
            if tend==False:
                tend = self.aim.data.t_all[-3]
            n = int(np.round((tend-t0)/dt))+1
            t_plot = np.linspace(t0,tend,n)
        
        elif calc==False:
            t_plot = self.aim.data.t_all

        return t_plot
    
    def clear_functions(self):
        try:
            del self.func
        except AttributeError:
            pass
        try:
            del self.sampled
        except AttributeError:
            pass
        
        return 0


    def make_functions(self,include=[],exclude=[],option='both',i='all',side=['l','r'],auto_clear=False,t=False,mode='mean_var',**kwargs):
        '''Obtains the returned properties'''
        import pointLISA
        from pointLISA import utils

        Nbins=kwargs.pop('Nbins',False)

        if auto_clear==True:
            self.clear_functions()

        ret=[]
        if include=='all' and exclude!='all':
            for k in OUTPUT.__dict__.keys():
                if 'get_' in k:
                    add = k.split('get_')[-1]
                    if add not in exclude:
                        ret.append(add)
        
        else:
            for inc in include:
                ret.append(inc)
         
        func=pointLISA.utils.Object() 
        sampled=pointLISA.utils.Object()
        
        if type(t)==bool:
            if t==False:
                t_plot = self.t_calc(calc=True) #add parameters
        else:
            t_plot=np.array(t)
        

        if i =='all':
            i=range(1,4)
        elif type(i)!=list:
            i=[i]
        if type(side)!=list:
            side=[side]

        for s in side:
            setattr(sampled,s,utils.Object())
            for i_sel in i:
                setattr(getattr(sampled,s),'i'+str(i_sel),utils.Object())



        ret_new = []
        for k in ret:
            ret_new.append(k.replace("'",""))
        ret = ret_new
        del ret_new

        for k in ret:
            print(k)
            
            if 'adjust' == k:
                for s in side:
                    for i_sel in i:
                        if self.aim.PAAM_deg==1:
                            if s=='l':
                                A = [np.array(getattr(self.aim,'t_adjust')[0][int(i_sel)-1])]
                                A.append(np.array(getattr(self.aim,'tele_adjust')[0][int(i_sel)-1])) 
                                #A.append(np.array([self.aim.tele_l_ang(i_sel,t) for t in A[0]]))
                            elif s=='r':
                                A = [np.array(getattr(self.aim,'t_adjust')[1][int(i_sel)-1])]
                                A.append(np.array(getattr(self.aim,'tele_adjust')[1][int(i_sel)-1])) 
                                #A.append(np.array([self.aim.tele_r_ang(i_sel,t) for t in A[0]]))
                        elif self.aim.PAAM_deg==2:
                            if s=='l':
                                A = [np.array(getattr(self.aim,'t_adjust')[0][int(i_sel)-1])]
                                A.append(np.array(getattr(self.aim,'tele_adjust')[0][int(i_sel)-1])) 
                            elif s=='r':
                                A = [np.array(getattr(self.aim,'t_adjust')[1][int(i_sel)-1])]
                                A.append(np.array(getattr(self.aim,'tele_adjust')[1][int(i_sel)-1]))   

#                            try:
#                                A = [np.array(getattr(self.aim,'t_adjust')[s][int(i_sel)])]
#                            except KeyError:
#                                A = [np.array(getattr(self.aim,'t_adjust')[str(i_sel)][s])]
#                            try:
#                                A.append(np.array(getattr(self.aim,'tele_adjust')[s][int(i_sel)]))
#                            except AttributeError:
#                                A.append(np.array(getattr(self.aim,'tele_ang_adjust')[str(i_sel)][s]))
                            
                        B = [A,'value='+k+', mode='+str(mode)]
                        setattr(getattr(getattr(sampled,s),'i'+str(i_sel)),k,B)
            else:
                if option=='both' or option=='function':
                    setattr(func,k,lambda i,t,side: self.mean_var(i,t,side,[k],Nbins=Nbins,mode=mode))
                if option=='both' or option=='sampled':
                    for s in side:
                        for i_sel in i:
                            A=[t_plot]
                            try:
                                A.append(np.array([self.mean_var(i_sel,t,s,ret=[k],Nbins=Nbins,mode=mode) for t in t_plot]))
                            except TypeError,e:
                                if "'dict' object is not callable" in str(e):
                                    A = getattr(self.aim,k)[s][i_sel]
                            B = [A,'value='+k+', mode='+str(mode)]
                            setattr(getattr(getattr(sampled,s),'i'+str(i_sel)),k,B)
        return [func,sampled]


### Calculate values/properties

def get_coor_tele(aim,i,t,side,tele_angle=False):
    '''Gets telescope coordinate system'''
    if tele_angle==False:
        if side == 'l':
            try:
                ret = aim.tele_l_coor(i,t)
            except AttributeError:
                tele_angle = aim.tele_l_ang(i,t)
        elif side =='r':
            try:
                ret = aim.tele_r_coor(i,t)
            except AttributeError:
                tele_angle = aim.tele_r_ang(i,t)

    try:
        return ret
    except:
        ret = methods.coor_tele(aim.data,i,t,tele_angle)
        return ret

def get_coor_beam_in__sun(aim,i,t,tdel,side,tele_angle_send=False,beam_angle_send=False,tele_angle_rec=False,offset=False,out=3):
    '''Gets incoming (received) beam coordinate system'''
    [i_self,i_left,i_right] = utils.i_slr(i)
    check=False
    if aim.data.calc_method=='Abram':
        tdel0 = 0
    elif aim.data.calc_method=='Waluschka':
        tdel0 = tdel
    try:
        if tele_angle_send==False and beam_angle_send==False:
            if side=='l':
                u_sun = aim.beam_r_coor(i_left,t-tdel)
            elif side=='r':
                u_sun = aim.beam_l_coor(i_right,t-tdel)
            check=True
    except AttributeError:
        check=False
        pass

    if check==False:
        if tele_angle_send is False:
            if side=='l':
                tele_angle_send = np.radians(30.0)
            elif side=='r':
                tele_angle_send = np.radians(-30.0)
        if beam_angle_send is False:
            if side=='l':
                beam_angle_send = 0.0
            elif side=='r':
                beam_angle_send = 0.0
        if offset is False:
            if side=='l':
                offset = get_offset(aim,i_left,t-tdel,'r')
            elif side=='r':
                offset = get_offset(aim,i_right,t-tdel,'l')
        elif offset == None:
            offset = 0.0

        if side=='l':
            u_sun = methods.beam_coor_out(aim.data,i_left,t-tdel,tele_angle_send,beam_angle_send,offset)
        elif side=='r':
            u_sun = methods.beam_coor_out(aim.data,i_right,t-tdel,tele_angle_send,beam_angle_send,offset)
    

    return u_sun

def get_coor_beam_out__send(aim,i,t,side,tele_angle=False,beam_angle=False,offset=False):
    '''Gets outgoing (transmitted) beam coordinate system'''
    check=False
    
    if check==False:
        if tele_angle is False:
            if side == 'l':
                tele_angle = aim.tele_l_ang(i,t)
            elif side =='r':
                tele_angle = aim.tele_r_ang(i,t)
        elif tele_angle ==None:
            if side == 'l':
                tele_angle = np.radians(-30.0)
            elif side =='r':
                tele_angle = np.radians(30.0)

        if beam_angle is False:
            if side == 'l':
                beam_angle = aim.beam_l_ang(i,t)
            elif side =='r':
                beam_angle = aim.beam_r_ang(i,t)
        elif beam_angle==None:
            if side == 'l':
                beam_angle = 0.0
            elif side =='r':
                beam_angle = 0.0
        
        if offset is False:
            offset = get_offset(aim,i,t,side)
        elif offset==None:
            offset=0.0

        ret = methods.beam_coor_out__send(aim.data,i,t,tele_angle,beam_angle,offset)
    return ret

def get_offset(aim,i,t,side):
    '''Gets the offset angle (inplane) between the telescope and transmitted beam inplane angle'''
    try:
        ret = aim.offset[side][i](t)
    except TypeError:
        try:
            ret = aim.offset[side][i] 
        except TypeError:
            try:
                ret = aim.offset(i,t,side)
            except TypeError:
                print(i,t,side)
                raise TypeError
    return ret

def get_start_calc(aim,i,t,side,tele_angle):
    '''Gets the starting point of the telescope (where the tranmitted beam is leaving the telescope'''
    try:
        if side=='l':
            ret = aim.tele_l_start(i,t)
        elif side=='r':
            ret = aim.tele_r_start(i,t)
    except AttributeError:
        ret = np.array(aim.data.putp(i,t)) + LA.unit(get_coor_tele(aim,i,t,side,tele_angle=tele_angle)[0])*aim.data.L_tele

    return ret

def values(inp,i,t,side,ksi=[0,0],mode='send',tele_angle_l=False,tele_angle_r=False,beam_angle_l=False,beam_angle_r=False,offset_l=False,offset_r=False,ret=[]):
    '''Runner function to obtain the output values for spacecraft i at time t'''
    from pointLISA import *
    [i_self,i_left,i_right] = utils.i_slr(i)
    

    if 'AIM' in str(inp):
        aim = inp
        outp = OUTPUT(aim)
    elif 'OUTPUT' in str(inp):
        aim = inp.aim
        outp = inp
    else:
        raise ValueError
     
    if tele_angle_l==None:
        tele_angle_l = -np.radians(30.0)
    if tele_angle_r==None:
        tele_angle_r = -np.radians(30.0)
    if beam_angle_l==None:
        beam_angle_l = 0.0
    if beam_angle_r==None:
        beam_angle_r = 0.0

    if mode=='send':
        if side=='l':
            tdel = aim.data.L_sl_func_tot(i_self,t)
            if aim.data.calc_method=='Waluschka':
                tdel0=tdel
            elif aim.data.calc_method=='Abram':
                tdel0=0            
            if offset_l is False:
                offset_l = get_offset(aim,i_self,t+tdel0,'l')
            elif offset_l == None:
                offset_l = 0.0
            if offset_r is False:
                offset_r = get_offset(aim,i_left,t+tdel,'r')
            elif offset_r == None:
                offset_r = 0.0
            i_send = i_self
            i_rec = i_left

        elif side=='r':
            tdel = aim.data.L_sr_func_tot(i_self,t)
            if aim.data.calc_method=='Waluschka':
                tdel0=tdel
            elif aim.data.calc_method=='Abram':
                tdel0=0
            if offset_l is False:
                offset_l = get_offset(aim,i_right,t+tdel,'l')
            elif offset_l == None:
                offset_l = 0.0
            if offset_r is False:
                offset_r = get_offset(aim,i_self,t+tdel0,'r')
            elif offset_r == None:
                offset_r = 0.0
            i_send = i_self
            i_rec = i_right

    elif mode=='rec':
        if side=='l':
            tdel = aim.data.L_rl_func_tot(i_self,t)
            if aim.data.calc_method=='Waluschka':
                tdel0=tdel
            elif aim.data.calc_method=='Abram':
                tdel0=0
            if offset_l is False:
                offset_l = get_offset(aim,i_self,t-tdel0,'l')
            elif offset_l==None:
                offset_l=0.0
            if offset_r is False:
                offset_r = get_offset(aim,i_left,t-tdel,'r')
            elif offset_r==None:
                offset_r=0.0
            i_send = i_left
            i_rec = i_self

        elif side=='r':
            tdel = aim.data.L_rr_func_tot(i_self,t)
            if aim.data.calc_method=='Waluschka':
                tdel0=tdel
            elif aim.data.calc_method=='Abram':
                tdel0=0
            if offset_l is False:
                offset_l = get_offset(aim,i_right,t-tdel,'l')
            elif offset_l==None:
                offset_l=0.0
            if offset_r is False:
                offset_r = get_offset(aim,i_self,t-tdel0,'r')
            elif offset_r==None:
                offset_r=0.0
            i_send = i_right
            i_rec = i_self

            
    if (mode=='send' and side=='l') or (mode=='rec' and side=='r'):
        tele_angle_start = tele_angle_l
        beam_angle_start = beam_angle_l
        tele_angle_end = tele_angle_r
        beam_angle_end = beam_angle_r
        offset_start = offset_l
        offset_end = offset_r
        
    elif (mode=='send' and side=='r') or (mode=='rec' and side=='l'):
        tele_angle_start = tele_angle_r
        beam_angle_start = beam_angle_r
        tele_angle_end = tele_angle_l
        beam_angle_end = beam_angle_l
        offset_start = offset_r
        offset_end = offset_l
 
#    if mode=='send': #...aberz
#        r_beam_start = methods.beam_coor_out(aim.data,i_send,t+tdel0,tele_angle_start,beam_angle_start,offset_start)[0]
#        tele_coor_start = methods.coor_tele(aim.data,i_send,t+tdel0,tele_angle_start)
#        r_beam_start_SC = methods.aberration_beam_coor(aim.data,i_send,t+tdel0,r_beam_start,reverse=True)
#        beam_angle_start = LA.angle(r_beam_start_SC,tele_coor_start[0])
#
#        r_beam_end = methods.beam_coor_out(aim.data,i_rec,t+tdel,tele_angle_end,beam_angle_end,offset_end)[0]
#        tele_coor_end = methods.coor_tele(aim.data,i_rec,t+tdel,tele_angle_end)
#        r_beam_end_SC = methods.aberration_beam_coor(aim.data,i_send,t+tdel,r_beam_start,reverse=True)
#        beam_angle_end = LA.angle(r_beam_end_SC,tele_coor_end[0])
#    
#    elif mode=='rec':
#        r_beam_start = methods.beam_coor_out(aim.data,i_send,t-tdel,tele_angle_start,beam_angle_start,offset_start)[0]
#        tele_coor_start = methods.coor_tele(aim.data,i_send,t-tdel,tele_angle_start)
#        r_beam_start_SC = methods.aberration_beam_coor(aim.data,i_send,t-tdel,r_beam_start,reverse=True)
#        beam_angle_start = LA.angle(r_beam_start_SC,tele_coor_start[0])
#
#        r_beam_end = methods.beam_coor_out(aim.data,i_rec,t-tdel0,tele_angle_end,beam_angle_end,offset_end)[0]
#        tele_coor_end = methods.coor_tele(aim.data,i_rec,t-tdel0,tele_angle_end)
#        r_beam_end_SC = methods.aberration_beam_coor(aim.data,i_send,t-tdel0,r_beam_end,reverse=True)
#        beam_angle_end = LA.angle(r_beam_end_SC,tele_coor_end[0])

    positions=utils.Object()
    positions.method = aim.data.calc_method
    positions.aim=aim
    positions.tele_angle_l = tele_angle_l
    positions.tele_angle_r = tele_angle_r
    positions.beam_angle_l = beam_angle_l
    positions.beam_angle_r = beam_angle_r
    positions.offset_l = offset_l
    positions.offset_r = offset_r
    
    param = ['mode','side','i_self','i_left','i_right','t','ksi','tdel','tdel0','tele_angle_start','tele_angle_end','beam_angle_start','beam_angle_end','aim','offset_start','offset_end','i_rec','i_send']
    for p in param:
        try:
            setattr(positions,p,locals()[p])
        except Exception as e:
            print(e)
            pass

    for r in ret:
        if r not in positions.__dict__.keys():
            try:
                positions_new = getattr(outp,'get_'+r)(positions)
                del positions
                positions = positions_new
            except AttributeError,e:
                print(e) 
                setattr(positions,r,getattr(aim,r)(i,t))
    
    return positions

def aberration(pos,u,mode='rec',**kwargs): # Only classical  
    if 'Object' in str(type(pos)):
        aber = pos.aim.data.aberration
    elif 'STATIC' in str(pos):
        aber = pos.data.aberration

    if aber==True:
        if 'Object' in str(type(pos)):
            if pos.aim.data.aberration==False:
                u_new = u
            elif pos.aim.data.aberration==True:
                if mode=='rec':
                    i = pos.i_rec
                    if pos.mode=='rec':
                        t = pos.t-pos.tdel0
                    elif pos.mode=='send':
                        t = pos.t-pos.tdel
                elif mode=='send':
                    i = pos.i_send
                    if pos.mode=='rec':
                        t = pos.t+pos.tdel
                    elif pos.mode=='send':
                        t = pos.t+pos.tdel0
        elif 'STATIC' in str(pos):
            t = kwargs['t']
            i = kwargs['i']

        V = -pos.aim.data.vel.abs(i,t)
        u_mag = np.linalg.norm(u)
        c_vec = LA.unit(u)*c
        u_new = LA.unit(c_vec+V)*u_mag
    
    elif aber==False:
        u_new = u
    
    return u_new


# Calculating the PAAM and telescope poiting angles
def tele_center_calc(aim,i,t,scale=1,lim=1e-12,max_count=5,print_on=False,value=0,tele_l=False,tele_r=False,beam_l=False,beam_r=False,offset_l=False,offset_r=False):
    '''Obtains the telescope pointing angle when the telesope is pointed with the center method'''
    [i_self,i_left,i_right] = utils.i_slr(i)
    
    lim = np.radians(5.0)
    if tele_l is False:
        tele_l=aim.tele_l_ang(i_self,t)
    elif tele_l==None:
        tele_l=np.radians(np.float64(-30.0))
    if tele_r is False:
        tele_r=aim.tele_r_ang(i_left,t)
    elif tele_r==None:
        tele_r=np.radians(np.float64(30.0))
    if beam_l is False:
        beam_l=aim.beam_l_ang(i_self,t)
    elif beam_l==None:
        beam_l=np.float64(0.0)
    if beam_r is False:
        beam_r=aim.beam_r_ang(i_self,t)
    elif beam_r==None:
        beam_r=np.float64(0.0)
    
    tele_l_old = tele_l
    tele_r_old = tele_r
    
    pos_send = lambda tele_l: values(aim,i_self,t,'l',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['off']).xoff
    send_solve = lambda tele_l: pos_send(tele_l)-value
    

    try:
        tele_l_new = scipy.optimize.brentq(send_solve,-lim-np.radians(30.0),lim-np.radians(30.0))
    except ValueError,e:
        if str(e)=='f(a) and f(b) must have different signs':
            tele_l_new=np.nan
 
    if tele_l_new!=np.nan:
        pos_rec = lambda tele_r: values(aim,i_left,t,'r',tele_angle_l=tele_l_new,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['off']).xoff
        rec_solve = lambda tele_r: pos_rec(tele_r)-value
        
        try:
            tele_r_new = scipy.optimize.brentq(rec_solve,-lim+np.radians(30.0),lim+np.radians(30.0))
        except ValueError,e:
            if str(e)=='f(a) and f(b) must have different signs':
                tele_r_new=np.nan
    else:
        tele_r_new = np.nan
        
    return [[tele_l_new,tele_r_new], False]

def PAAM_center_calc(aim,i,t,para='yoff_ab',scale=1,lim=1e-12,max_count=5,print_on=False,tele_l=None,tele_r=None,beam_l=None,beam_r=None,offset_l=None,offset_r=None,value=0,method='iter',margin=0.01):
    '''Obtains the PAAM pointing angle when the PAAM is pointed with the center method'''
    [i_self,i_left,i_right] = utils.i_slr(i)
    
    if tele_l is False:
        tele_l = aim.tele_l_ang(i,t)
    elif tele_l==None:
        tele_l = np.radians(np.float64(-30.0))
    if tele_r is False:
        tele_r = aim.tele_r_ang(i,t)
    elif tele_r == None:
        tele_r = np.radians(np.float64(30.0))
    if beam_l is False:
        beam_l=aim.beam_l_ang(i,t)
    elif beam_l==None:
        beam_l = np.float64(0.0)
    if beam_r is False:
        beam_r=aim.beam_r_ang(i,t)
    elif beam_r==None:
        beam_r=np.float64(0.0)
    
    lim = 1.0e-3
    beam_l_old = beam_l
    beam_r_old = beam_r
    pos_send = lambda beam_l: values(aim,i_self,t,'l',tele_angle_l=tele_l,tele_angle_r=tele_l,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['off']).yoff
    send_solve = lambda beam_l: pos_send(beam_l)-value


    try:
        beam_l_new = scipy.optimize.brentq(send_solve,-lim,lim)
    except ValueError,e:
        if str(e)=='f(a) and f(b) must have different signs':
            beam_l_new=np.nan

    if beam_l_new!=np.nan:
        pos_rec = lambda beam_r: values(aim,i_left,t,'r',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam_l_new,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['off']).yoff
        rec_solve = lambda beam_r: pos_rec(beam_r)-value

        try:
            beam_r_new = scipy.optimize.brentq(rec_solve,-lim,lim)
        except ValueError,e:
            if str(e)=='f(a) and f(b) must have different signs':
                beam_r_new=np.nan
    else:
        beam_r_new = np.nan

    mode='Converged'

    return [[beam_l_new,beam_r_new], mode]

def PAAM_wavefront_calc(aim,i,t,side,lim=1e-9,tele_l=False,tele_r=False):
    '''Obtains the PAAM pointing angle when the PAAM is pointed with the wavefront method'''
    if side=='l':
        angy = lambda beam: values(aim,i,t,'l',mode='send',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam,ret=['angy_wf_rec']).angy_wf_rec
    elif side=='r':
        angy = lambda beam: values(aim,i,t,'r',mode='send',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_r=beam,ret=['angy_wf_rec']).angy_wf_rec

    try:
        ret = scipy.optimize.brentq(angy,-1e-5,1e-5,xtol=lim)
    except ValueError,e:
        if str(e)=='f(a) and f(b) must have different signs':
            ret=np.nan
    return ret



def tele_wavefront_calc(aim,i_l,t,scale=1,lim=1e-12,max_count=5,print_on=False,value=0,tele_l=False,tele_r=False,beam_l=False,beam_r=False,offset_l=False,offset_r=False):
    '''Obtains the telescope pointing angle when the telesope is pointed with the center method'''
    [i_self,i_left,i_right] = utils.i_slr(i_l)

    lim = np.radians(5.0)
    if tele_l is False:
        tele_l=aim.tele_l_ang(i_self,t)
    elif tele_l==None:
        tele_l=np.radians(np.float64(-30.0))
    if tele_r is False:
        tele_r=aim.tele_r_ang(i_left,t)
    elif tele_r==None:
        tele_r=np.radians(np.float64(30.0))
    if beam_l is False:
        beam_l=aim.beam_l_ang(i_self,t)
    elif beam_l==None:
        beam_l=np.float64(0.0)
    if beam_r is False:
        beam_r=aim.beam_r_ang(i_self,t)
    elif beam_r==None:
        beam_r=np.float64(0.0)
    
    tele_l_old = tele_l
    tele_r_old = tele_r
    
    par = 'angx_wf_rec'

    pos_l = getattr(values(aim,i_self,t,'l',mode='rec',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=[par]),par)
    tele_l_new = tele_l-pos_l
    
    pos_r = getattr(values(aim,i_left,t,'r',mode='rec',tele_angle_l=tele_l_new,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=[par]),par)
    tele_r_new = tele_r-pos_r
    
    return [[tele_l_new,tele_r_new], [pos_l,pos_r]]


#def tele_wavefront_calc(aim,i_l,t,method,scale=1,lim=1e-12,max_count=20,tele_angle_l=None,tele_angle_r=None,beam_l=None,beam_r=None,offset_l=False,offset_r=False,print_on=False): #The sending is always the left telescope and the receiving the right one
#    '''Obtains the telescope pointing angle when the telesope is pointed with the wavefront method'''
#    import pointLISA.utils
#
#    i_r = pointLISA.utils.i_slr(i_l)[1]
#    
#    tdel = aim.data.L_sl_func_tot(i_l,t)
#    t_l=t
#    t_r=t+tdel #...because considering traveling time makes it more complex (solve for functions)
#
#    if tele_angle_l is False:
#        tele_angle_l=aim.tele_l_ang(i_l,t_l)
#    elif tele_angle_l == None:
#        tele_angle_l=np.radians(-30.0)
#    if tele_angle_r is False:
#        tele_angle_r=aim.tele_r_ang(i_r,t_r)
#    elif tele_angle_r ==None:
#        tele_angle_r = np.radians(30.0)
#    
#    calc_l=100
#    calc_r=100
#    calc_l_old=calc_l
#    calc_r_old=calc_r
#
#    count=0
#    lim_val=100
#
#    para = 'angx_wf_send'
#    
#    tdel = aim.data.L_sl_func_tot(i_l,t)
#
#    if method=='iter':
#        while count<max_count:
#            calc_l_old = calc_l
#            calc_r_old = calc_r
#            tele_angle_l_old = tele_angle_l
#            tele_angle_r_old = tele_angle_r
#            pos_l = values(aim,i_l,t_l,'l',mode='send',tele_angle_l=tele_angle_l,tele_angle_r=tele_angle_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r)
#            calc_l = getattr(getattr(OUTPUT(aim),'get_'+para)(pos_l),para)
#            tele_angle_r = tele_angle_r+scale*calc_l
#
#            pos_r = values(aim,i_r,t_r+tdel,'r',mode='send',tele_angle_l=tele_angle_l,tele_angle_r=tele_angle_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r)
#            calc_r = getattr(getattr(OUTPUT(aim),'get_'+para)(pos_r),para)
#            tele_angle_l = tele_angle_l+scale*calc_r
#            count=count+1
#            lim_val = max(abs(abs(calc_l_old)-abs(calc_l)),abs(abs(calc_r_old)-abs(calc_r)))
#            if print_on:
#                print(tele_angle_l,tele_angle_r)
#                print(calc_l,calc_r)
#                print(max(abs(calc_l),abs(calc_r)),max(abs(calc_l),abs(calc_r))<lim)
#                print(lim_val)
#                print('')
#
#            if lim_val<=lim:
#                mode = 'Result is converged'
#                if print_on:
#                    print(mode)
#                break
#
#            if count>= max_count:
#                mode = 'Maximum iteration limit has been reached'
#                if print_on:
#                    print(mode)
#                break
#    return [[tele_angle_l,tele_angle_r],mode]



def get_tele_wavefront(aim,i,t,side,method,scale=1,lim=1e-12,max_count=20,print_on=False,value=0.0,tele_angle_l=None,tele_angle_r=None,beam_l=None,beam_r=None,offset_l=False,offset_r=False): 
    '''Gets all telescope pointing angles along one arm for the wavefront method'''
    if side=='l':
        i_l = i
        tdel=0
    elif side=='r':
        i_l = utils.i_slr(i)[2]
        i_r = i
        tdel = aim.data.L_rr_func_tot(i_r,t)
    
    tele_l_old = 0.0
    tele_r_old = 0.0
    if tele_angle_l==None:
        tele_l = np.radians(-30.0)
    else:
        tele_l = tele_angle_l
    if tele_angle_r==None:
        tele_r = np.radians(30.0)
    else:
        tele_r = tele_angle_r

    
    count=0
    while count<max_count:
        [[tele_l_new,tele_r_new],con] = tele_wavefront_calc(aim,i,t,tele_l=tele_l,tele_r=tele_r,beam_l=beam_l,beam_r=beam_r,offset_l=offset_l,offset_r=offset_r)
        count = count+1
        if count>= max_count:
            mode = 'Maximum iteration limit has been reached'
            tele_l = tele_l_new
            tele_r = tele_r_new           
            if print_on:
                print(mode)
            break
        elif max(con)<1.0e-9: #max(tele_l_new-tele_l,tele_r_new-tele_r)<1.0e9:
            mode = 'Result is converged'
            tele_l = tele_l_new
            tele_r = tele_r_new
            if print_on:
                print(mode)
            break
        tele_l = tele_l_new
        tele_r = tele_r_new
        
    if side=='l':
        return tele_l

    elif side=='r':
        return tele_r



    
