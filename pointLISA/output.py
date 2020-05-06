from pointLISA import *

class OUTPUT():
    '''This class uses an AIM object and calculate some useful properties (output) as attributes of an OUTPUT object'''
    def __init__(self,aim=False,**kwargs):
        if aim!=False:
            self.aim = aim
        else:
            self.aim = utils.Object()
            self.aim.data = kwargs['data']

        #### Obtain parameters
        #for k, value in parameters.__dict__.items():
        #    if '__' not in k:
        #        globals()[k] = value
        #        setattr(self,k,value)
    
    def pupil(self,Nbins=2): #Aperture
        '''Creates pixels on the aperture which can be useful to examine the variation of a property over the aperture'''
        try:
            del_at = ['xlist','ylist','Deltax','Deltay','Nbinsx','Nbinsy']
            for at in del_at: 
                delattr(self,at)
        except:
            pass

        D_calc=self.aim.data.param.D

        xlist = np.linspace(-D_calc*0.5,D_calc*0.5,Nbins+1)
        self.xlist = xlist[0:-1]+0.5*(xlist[1]-xlist[0])
        self.ylist = self.xlist
        self.Deltax = self.xlist[1]-self.xlist[0]
        self.Deltay = self.ylist[1]-self.ylist[0]
        self.Nbinsx = len(self.xlist)
        self.Nbinsy = len(self.ylist)

    def w(self,z):
        '''Beamwaist as a function of z (z=coordinate along beamline)'''
        zR = np.pi*(self.aim.data.param.w0_laser**2)/self.aim.data.param.labda

        return self.aim.data.param.w0_laser*((1+((z/zR)**2))**0.5)

    def R(self,z,guess=False):
        '''The radius of curvasture R as a function of z'''
        if z!=np.nan:
            zR = np.pi*(pos.aim.data.param.w0_laser**2)/self.labda

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
                    ret=calc.get_coor_tele(pos.aim,pos.i_self,pos.t,pos.side,tele_angle=pos.tele_angle_start)
                elif pos.mode=='rec':
                    ret=calc.get_coor_tele(pos.aim,pos.i_opp,pos.t-pos.tdel,pos.invside,tele_angle=pos.tele_angle_start)
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
                    ret=calc.get_coor_beam_out__send(pos.aim,pos.i_send,pos.t,pos.side,tele_angle=pos.tele_angle_start,beam_angle=pos.beam_angle_start,offset=pos.offset_start)
                elif pos.mode=='rec':
                    ret=calc.get_coor_beam_out__send(pos.aim,pos.i_send,pos.t-pos.tdel,pos.invside,tele_angle=pos.tele_angle_start,beam_angle=pos.beam_angle_start,offset=pos.offset_start)
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
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t+pos.tdel0,v,reverse=False)
                elif pos.mode=='rec':
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t-pos.tdel,v,reverse=False)
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
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t+pos.tdel,v,reverse=True)
                elif pos.mode=='rec':
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t-pos.tdel0,v,reverse=True)
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
                    ret = calc.get_coor_beam_out__send(pos.aim,pos.i_send,pos.t+pos.tdel,pos.side,tele_angle=pos.tele_angle_start,beam_angle=pos.beam_angle_start,offset=pos.offset_start)
                elif pos.mode=='rec':
                    ret = calc.get_coor_beam_out__send(pos.aim,pos.i_send,pos.t-pos.tdel0,pos.invside,tele_angle=pos.tele_angle_start,beam_angle=pos.beam_angle_start,offset=pos.offset_start)
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
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t+pos.tdel0,v,reverse=False)
                elif pos.mode=='rec':
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t-pos.tdel,v,reverse=False)
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
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t+pos.tdel,v,reverse=True)
                elif pos.mode=='rec':
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t-pos.tdel0,v,reverse=True)
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
                    ret = calc.get_coor_tele(pos.aim,pos.i_opp,pos.t+pos.tdel,pos.invside,tele_angle=pos.tele_angle_end)
                elif pos.mode=='rec':
                    ret = calc.get_coor_tele(pos.aim,pos.i_self,pos.t-pos.tdel0,pos.side,tele_angle=pos.tele_angle_end)
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
                    ret = calc.get_start_calc(pos.aim,pos.i_self,pos.t+pos.tdel0,pos.side,pos.tele_angle_start)
                elif pos.mode=='rec':
                    ret = calc.get_start_calc(pos.aim,pos.i_opp,pos.t-pos.tdel,pos.invside,pos.tele_angle_start)
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
                    ret = calc.get_start_calc(pos.aim,pos.i_opp,pos.t+pos.tdel,pos.invside,pos.tele_angle_end)
                elif pos.mode=='rec':
                    ret = calc.get_start_calc(pos.aim,pos.i_self,pos.t-pos.tdel0,pos.side,pos.tele_angle_end)
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
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t+pos.tdel0,v,reverse=True)
                elif pos.mode=='rec':
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t-pos.tdel,v,reverse=True)            
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
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t+pos.tdel,v,reverse=True)
                elif pos.mode=='rec':
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t-pos.tdel0,v,reverse=True)            
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
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t+pos.tdel,pos.startend__sun,reverse=True)
                elif pos.mode=='rec':
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t-pos.tdel0,pos.startend__sun,reverse=True)
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
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t+pos.tdel0,pos.startend__sun,reverse=True) #.......check on time it emits!!!
                elif pos.mode=='rec':
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t-pos.tdel,pos.startend__sun,reverse=True)
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
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t+pos.tdel0,v,reverse=False)
                elif pos.mode=='rec':
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t-pos.tdel,v,reverse=False)
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
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t+pos.tdel,v,reverse=True)
                elif pos.mode=='rec':
                    ret = calc.aberration_beam_coor(pos.aim.data,pos.i_rec,pos.t-pos.tdel0,v,reverse=True)
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
                    v_send = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t+pos.tdel0,v_sun,reverse=True)
                elif pos.mode=='rec':
                    v_send = calc.aberration_beam_coor(pos.aim.data,pos.i_send,pos.t-pos.tdel,v_sun,reverse=True)
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
                if pos.alpha>pos.aim.data.param.FOV:#pos.waist>pos.aim.data.FOV:
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
                ret = (pos.I0*np.exp((-2*(pos.xoff**2+pos.yoff**2))/(pos.waist**2)))*(pos.aim.data.param.w0_laser/pos.waist)

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
                ret = (pos.I0*np.exp((-2*(pos.xoff**2))/(pos.waist**2)))*(pos.aim.data.param.w0_laser/pos.waist)

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
                ret = (pos.aim.data.param.P_L*np.pi*(pos.aim.data.param.w0_laser**2))/2.0
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
#                        ret = calc.get_coor_tele(pos.aim,pos.i_self,pos.t+pos.tdel0,'l',tele_angle=pos.tele_angle_l)
#                    elif pos.side=='r':
#                        ret = calc.get_coor_tele(pos.aim,pos.i_self,pos.t+pos.tdel0,'r',tele_angle=pos.tele_angle_r)
#                elif pos.mode=='rec':
#                    if pos.side=='l':
#                        ret = calc.get_coor_tele(pos.aim,pos.i_left,pos.t-pos.tdel,'r',tele_angle=pos.tele_angle_r)
#                    elif pos.side=='r':
#                        ret = calc.get_coor_tele(pos.aim,pos.i_right,pos.t-pos.tdel,'l',tele_angle=pos.tele_angle_l)
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
                        ret = calc.get_coor_tele(pos.aim,pos.i_left,pos.t+pos.tdel,'r',tele_angle=pos.tele_angle_r)
                    elif pos.side=='r':
                        ret = calc.get_coor_tele(pos.aim,pos.i_right,pos.t+pos.tdel,'l',tele_angle=pos.tele_angle_l)
                elif pos.mode=='rec':
                    if pos.side=='l':
                        ret = calc.get_coor_tele(pos.aim,pos.i_self,pos.t-pos.tdel0,'l',tele_angle=pos.tele_angle_l)
                    elif pos.side=='r':
                        ret = calc.get_coor_tele(pos.aim,pos.i_self,pos.t-pos.tdel0,'r',tele_angle=pos.tele_angle_r)
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
#                        ret = calc.get_start_calc(pos.aim,pos.i_self,pos.t+pos.tdel0,'l',pos.tele_angle_l)
#                    elif pos.side=='r':
#                        ret = calc.get_start_calc(pos.aim,pos.i_self,pos.t+pos.tdel0,'r',pos.tele_angle_r)
#                elif pos.mode=='rec':
#                    if pos.side=='l':
#                        ret = calc.get_start_calc(pos.aim,pos.i_left,pos.t-pos.tdel,'r',pos.tele_angle_r)
#                    elif pos.side=='r':
#                        ret = calc.get_start_calc(pos.aim,pos.i_right,pos.t-pos.tdel,'l',pos.tele_angle_l)
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
#                        ret = calc.get_start_calc(pos.aim,pos.i_left,pos.t+pos.tdel,'r',pos.tele_angle_r)+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
#                    elif pos.side=='r':
#                        ret = calc.get_start_calc(pos.aim,pos.i_right,pos.t+pos.tdel,'l',pos.tele_angle_l)+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
#                elif pos.mode=='rec':
#                    if pos.side=='l':
#                        ret = calc.get_start_calc(pos.aim,pos.i_self,pos.t-pos.tdel0,'l',pos.tele_angle_l)+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
#                    elif pos.side=='r':
#                        ret = calc.get_start_calc(pos.aim,pos.i_self,pos.t-pos.tdel0,'r',pos.tele_angle_r)+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
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

        func = lambda x,y: calc.values(self,i,t,side,ret=ret,ksi=[x,y],tele_angle_l=tele_angle_l,tele_angle_r=tele_angle_r,beam_angle_l=beam_angle_l,beam_angle_r=beam_angle_r)

        if mode=='center':
            return getattr(func(0,0),ret[0])
        elif 'var' in mode or 'mean' in mode:
            A=[]
            for x in self.xlist:
                for y in self.ylist:
                    if ((x**2+y**2)**0.5)>=self.aim.data.param.D/2.0:
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

#OUTPUT = calculations_output

