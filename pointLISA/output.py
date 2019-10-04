from imports import *
import inspect
import numpy as np
import methods

class OUTPUT():
    def __init__(self,aim=False,**kwargs):
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
        try:
            del_at = ['xlist','ylist','Deltax','Deltay','Nbinsx','Nbinsy']
            for at in del_at: 
                delattr(self,at)
            #del self.xlist,self.ylist,self,Deltax,self.Deltay,self.Nbinsx,self.Nbinsy
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

    def w(self,z): #Beamwaist as a function of z
        zR = np.pi*(self.w0_laser**2)/self.labda

        return self.w0_laser*((1+((z/zR)**2))**0.5)

    def R(self,z,guess=False):
        if z!=np.nan:
            zR = np.pi*(self.w0_laser**2)/self.labda

            if guess==False:
                return abs(z*(1+((zR/z)**2)))

            elif guess==True:
                return z
        else:
            return np.nan
    
    def w0(self,angx,angy,ksi):
        # Opposite wfe
        #if side=='l':
        #    angy = self.aim.beam_l_ang(i,t)
        #elif side == 'r':
        #    angy = self.aim.beam_r_ang(i,t)
        #angx=0#...add jitter

        [x,y] = ksi

        labda = self.labda
        w_error = 2*np.pi*((x*np.sin(angx)+y*np.sin(angy))/labda)

        return w_error
    
    def u0(self,ksi):#...for gaussian beam, adjust when implementing other beams
        w = self.w(0.0)
        [x,y] = ksi

        return np.exp(-((x**2+y**2)/(w**2)))



    def aperture(self,xlist,ylist,function,dType=np.float64): # Creates matrix of function over an aperture (circle)
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
        try:
            if z!=np.nan:
                x = np.float64(x)
                y = np.float64(y)
                z = np.float64(z)
                #lim = (x**2+y**2) 
                #R_new = lambda dz: self.R(z+dz,guess=False)
                f_solve = lambda dz: (self.R(z+dz,guess=R_guess) - (self.R(z+dz,guess=R_guess)**2 - (x**2+y**2))**0.5) - dz
                f_solve_2 = lambda dz: (z- (((z+dz)**2 - x**2 -y**2 )**0.5))
                dz_sol = scipy.optimize.brentq(f_solve,-0.5*z,0.5*z,xtol=1e-64)
                #dz_sol_3 = scipy.optimize.brentq(f_solve_2,-lim,lim,xtol=1e-64)
                #dz_sol_2 = scipy.optimize.fsolve(f_solve,0,xtol=1e-128)
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
                    #print(dz_sol,1e-64)
                    #print(x,y)
                    #print(dz_sol,dz_sol_2,dz_sol_3)
                    #print(x,y,z)
                    return [z+dz_sol,dz_sol]

        except RuntimeError:
            #print(x,y,z)
            if ret=='piston':
                return np.nan
            elif ret=='all':
                return [np.nan,np.nan]


    def get_coordinate_systems(self,speed=0):
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
        func = 'get_'+str(e).split("attribute '")[-1].replace("'","")
        pos = getattr(self,func)(pos)

        return pos

    ### Returns important parameters
     
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

    def get_offset(self,pos):
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

    def get_off(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_startbeam_out,pos.end-pos.start)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                setattr(pos,'xoff',ret[2])
                setattr(pos,'yoff',ret[1])
                setattr(pos,'zoff',ret[0])
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos
    
    def get_tele_ang(self,pos):
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

    def get_PAAM_ang(self,pos):
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
    
### begin 16 september
    def get_invside(self,pos):
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

    def get_i_opp(self,pos):
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

    def get_coor_startbeam(self,pos):
        check=False
        while check==False:
            try:
                if pos.mode=='send':
                    ret=get_coor_beam_out(pos.aim,pos.i_self,pos.t,pos.side,tele_angle=pos.tele_angle_start,beam_angle=pos.beam_angle_start,offset=pos.offset_start)
                elif pos.mode=='rec':
                    ret=get_coor_tele(pos.aim,pos.i_opp,pos.t-pos.tdel,pos.invside,tele_angle=pos.tele_angle_start,beam_angle=pos.beam_angle_start)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_vecs_endbeam(self,pos):
        check=False
        while check==False:
            try:
                if pos.mode=='send':
                    ret=get_coor_beam_in(pos.aim,pos.i_opp,pos.t+pos.tdel,pos.tdel,pos.invside,tele_angle_send=pos.tele_angle_start,beam_angle_send=pos.beam_angle_start,tele_angle_rec=pos.tele_angle_end,offset=pos.offset_start)

                elif pos.mode=='rec':
                    ret=get_coor_beam_in(pos.aim,pos.i_self,pos.t,pos.tdel,pos.side,tele_angle_send=pos.tele_angle_start,beam_angle_send=pos.beam_angle_start,tele_angle_rec=pos.tele_angle_end,offset=pos.offset_start)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_vec_endbeam(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.vecs_endbeam[0]
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_vec_endbeam_na(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.vecs_endbeam[1]
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_vec_endbeam_ab(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.vecs_endbeam[2]
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
                    ret = get_coor_tele(pos.aim,pos.i_opp,pos.t+pos.tdel,pos.invside,tele_angle=tele_angle_end)
                elif pos.mode=='rec':
                    ret = get_coor_tele(pos.aim,pos.i_self,pos.t,pos.side,tele_angle=tele_angle_end)
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
                    ret = get_start_calc(pos.aim,pos.i_awld,pos.t-pos.tdel0,pos.side,pos.tele_angle_end)
                if ksi!=[0,0]:
                    ret = ret+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos



### end 16 september

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


    def get_arm(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.end-pos.start
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos
    
    def get_arm_tele_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_end,-pos.arm)
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
                ret = LA.matmul(pos.coor_starttele,pos.arm)
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
                    #try:
                    #   [piston,z_extra] = self.z_solve(pos.xoff,pos.yoff,pos.zoff,ret='all')
                    #except:
                    #    [piston,z_extra] = [np.nan,np.nan]
                    ret = self.R(pos.piston)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos) 
        return pos

    def get_R_vec_beam_send(self,pos):
        check=False
        while check==False:
            try:
                ret = np.array([(pos.R**2-pos.xoff**2-pos.yoff**2)**0.5,pos.yoff,pos.xoff]) # In beam frame
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
                ret = LA.matmul(pos.coor_starttele,LA.matmul(np.linalg.inv(pos.coor_startbeam_out),pos.R_vec_beam_send))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_R_vec_origin(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(np.linalg.inv(pos.coor_startbeam_out),pos.R_vec_beam_send)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_R_vec_tele_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_end,pos.R_vec_origin)
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
                ret = LA.matmul(pos.coor_end,-pos.vec_endbeam)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos
    
    def get_beam_receive_rec_na(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_end,-pos.vec_endbeam_na)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angx_ab_rec(self,pos):
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

    def get_angy_ab_rec(self,pos):
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

    def get_angx_nab_rec(self,pos):
        check=False
        while check==False:
            try:

                ret = np.sign(pos.beam_receive_rec_na[2])*abs(np.arctan(pos.beam_receive_rec_na[2]/pos.beam_receive_rec_na[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angy_nab_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.beam_receive_rec_na[1])*abs(np.arctan(pos.beam_receive_rec_na[1]/pos.beam_receive_rec_na[0]))
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
                ret = LA.matmul(pos.coor_starttele,pos.vec_endbeam)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos
    
    def get_beam_receive_send_na(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_starttele,pos.vec_endbeam_na)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angx_ab_send(self,pos):
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

    def get_angy_ab_send(self,pos):
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

    def get_angx_nab_send(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.beam_receive_send_na[2])*abs(np.arctan(pos.beam_receive_send_na[2]/pos.beam_receive_send_na[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angy_nab_send(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.beam_receive_send_na[1])*abs(np.arctan(pos.beam_receive_send_na[1]/pos.beam_receive_send_na[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

###


#    def get_tele_ab_rec(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = LA.matmul(pos.coor_starttele,LA.matmul(np.linalg.inverse(pos.coor_end),pos.beam_receive_rec))
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos
#
#    def get_tele_nab_rec(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = LA.matmul(pos.coor_starttele,LA.matmul(np.linalg.inverse(pos.coor_end),pos.beam_receive_rec_na))
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos
#
#    def get_angx_ab_send(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = np.sign(pos.tele_ab_rec[2])*abs(np.arctan(pos.tele_ab_rec[2]/pos.tele_ab_rec[0]))
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos
# 
#    def get_angy_ab_send(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = np.sign(pos.tele_ab_rec[1])*abs(np.arctan(pos.tele_ab_rec[1]/pos.tele_ab_rec[0]))
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos
#
#    def get_angx_nab_send(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = np.sign(pos.tele_nab_rec[2])*abs(np.arctan(pos.tele_nab_rec[2]/pos.tele_nab_rec[0]))
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos
# 
#    def get_angy_nab_send(self,pos):
#        check=False
#        while check==False:
#            try:
#                ret = np.sign(pos.tele_nab_rec[1])*abs(np.arctan(pos.tele_nab_rec[1]/pos.tele_nab_rec[0]))
#                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
#                check=True
#            except AttributeError,e:
#                #print(e)
#                self.add_attribute(e,pos)
#        return pos

    def get_R_vec_tele_rec(self,pos): #This vector is reversed (pointed away from receiving telescope)
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_end,-pos.R_vec_origin)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos


  
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

    def get_angx_R(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.R_vec_beam_send[2])*abs(np.arctan(pos.R_vec_beam_send[2]/pos.R_vec_beam_send[0]))
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
                ret = np.sign(pos.R_vec_beam_send[1])*abs(np.arctan(pos.R_vec_beam_send[1]/pos.R_vec_beam_send[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angx_tele(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.tele_vec[2])*abs(np.arctan(pos.tele_vec[2]/pos.tele_vec[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angy_tele(self,pos):
        check=False
        while check==False:
            try:
                ret = np.sign(pos.tele_vec[1])*abs(np.arctan(pos.tele_vec[1]/pos.tele_vec[0]))
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
                #ret = (pos.angx_tele-pos.angx_R)
                ret = np.arctan(abs(pos.R_vec_tele_rec[2]/pos.R_vec_tele_rec[0]))*np.sign(pos.R_vec_tele_rec[2])+pos.angx_ab_rec-pos.angx_nab_rec
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
                ret = np.arctan(abs(pos.R_vec_tele_rec[1]/pos.R_vec_tele_rec[0]))*np.sign(pos.R_vec_tele_rec[1])+pos.angy_ab_rec-pos.angy_nab_rec
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
                ret = np.arctan(abs(pos.R_vec_tele_send[2]/pos.R_vec_tele_send[0]))*np.sign(pos.R_vec_tele_send[2])+pos.angx_ab_send-pos.angx_nab_send
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
                ret = np.arctan(abs(pos.R_vec_tele_send[1]/pos.R_vec_tele_send[0]))*np.sign(pos.R_vec_tele_send[1])+pos.angy_ab_send-pos.angy_nab_send
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
                tele_end = LA.matmul(pos.coor_starttele,-pos.coor_end[0]) # tele_rec is reversed
                angx_tele_rec = np.sign(tele_end[2])*np.arctan(abs(tele_end[2]/tele_end[0]))
                angx_R_vec_tele_send = np.sign(pos.R_vec_tele_send[2])*np.arctan(abs(pos.R_vec_tele_send[2]/pos.R_vec_tele_send[0]))
                ret = angx_tele_rec - angx_R_vec_tele_send+pos.angx_ab_send-pos.angx_nab_send
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
                tele_end = LA.matmul(pos.coor_starttele,pos.coor_end[0])
                R_vec_tele_send_y = pos.R_vec_tele_send
                R_vec_tele_send_y[2] = 0
                R_y_origin = LA.matmul(np.linalg.inv(pos.coor_startbeam),R_vec_tele_send_y)
                R_y_rec = LA.matmul(pos.coor_end,-R_y_origin)
                ret = np.sign(R_y_rec[1])*np.arctan(abs(R_y_rec[1]/R_y_rec[0]))+pos.angy_ab_send-pos.angy_nab_send
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angx_wf_rec(self,pos):
        check=False
        while check==False:
            try:
                R_vec_rec = -pos.R_vec_tele_rec
                ret = np.sign(pos.R_vec_tele_rec[2])*np.arctan(abs(pos.R_vec_tele_rec[2]/pos.R_vec_tele_rec[0]))
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
                R_vec_rec = -pos.R_vec_tele_rec
                ret = np.sign(pos.R_vec_tele_rec[1])*np.arctan(abs(pos.R_vec_tele_rec[1]/pos.R_vec_tele_rec[0]))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos



### Hier gebleven
    def get_beam_direction_origin(self,pos):
        check=False
        while check==False:
            try:
                ret = np.array(pos.coor_startbeam[0])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_beam_direction_rec(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_end,pos.beam_direction_origin)
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

    def get_FOV_beamline(self,pos):
        check=False
        while check==False:
            try:
                ret = np.arccos(-pos.beam_direction_rec[0]/np.linalg.norm(pos.beam_direction_rec))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_alpha(self,pos):
        check=False
        while check==False:
            try:
                ret = (abs(pos.angx_wf_send)**2 +abs(pos.angy_wf_send)**2)**0.5
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_FOV_position(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.angle(pos.start-pos.end,pos.coor_end[0])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_tilt(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.FOV_wavefront
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

#    def get_u(self,pos,Nbins_inp=2):
#        check=False
#        while check==False:
#            try:
#                r = pos.r
#                z=pos.zoff
#                angx=pos.angx_send
#                angy=pos.angy_send
#
#                try:
#                    Nbins=self.Nbins
#                except:
#                    self.Nbins = Nbins_inp
#                try:
#                    xlist=self.xlist
#                    ylist=self.ylist
#                except AttributeError:
#                    self.pupil(Nbins=self.Nbins)
#                    xlist=self.xlist
#                    ylist=self.ylist
#                labda = self.aim.data.labda
#                k = (2*np.pi)/labda
#
#                if len(xlist)==1 and len(ylist)==1:
#                    dksi = (self.D**2)*(np.pi/4.0)
#                else:
#                    dksi = (xlist[1]-xlist[0])*(ylist[1]-ylist[0])
#
#                ret=0
#                for i in range(0,len(xlist)):
#                    for j in range(0,len(ylist)):
#                        ksi = np.array([xlist[i],ylist[j]])
#                        T1 = np.exp((1j*k*np.dot(r,ksi))/z)
#                        T2 = self.u0(ksi)
#                        T3 = np.exp(1j*self.w0(angx,angy,ksi))
#                        ret = ret+T1*T2*T3
#                ret = ret*dksi*(1j*k*np.exp(-(1j*k*(np.linalg.norm(r)**2))/(2*z))/(2*np.pi*z))
#
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
                #ret = I_0*np.cos(pos.angx_rec)*np.cos(pos.angy_rec)
                #ret = (abs(pos.u)**2)[0]#*np.cos(pos.angx_rec)*np.cos(pos.angy_rec)
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

    def get_coor_starttele(self,pos):
        check=False
        while check==False:
            try:
                if pos.mode=='send':
                    if pos.side=='l':
                        ret = get_coor_tele(pos.aim,pos.i_self,pos.t,'l',tele_angle=pos.tele_angle_l)
                    elif pos.side=='r':
                        ret = get_coor_tele(pos.aim,pos.i_self,pos.t,'r',tele_angle=pos.tele_angle_r)
                elif pos.mode=='rec':
                    if pos.side=='l':
                        ret = get_coor_tele(pos.aim,pos.i_left,pos.t-pos.tdel,'r',tele_angle=pos.tele_angle_r)
                    elif pos.side=='r':
                        ret = get_coor_tele(pos.aim,pos.i_right,pos.t-pos.tdel,'l',tele_angle=pos.tele_angle_l)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_coor_startbeam_out(self,pos):
        check=False
        while check==False:
            try:
                if pos.mode=='send':
                    if pos.side=='l':
                        ret = get_coor_beam_out(pos.aim,pos.i_self,pos.t,'l',tele_angle=pos.tele_angle_l,beam_angle=pos.beam_angle_l,offset=pos.offset_l)
                        
                    elif pos.side=='r':
                        ret = get_coor_beam_out(pos.aim,pos.i_self,pos.t,'r',tele_angle=pos.tele_angle_r,beam_angle=pos.beam_angle_r,offset=pos.offset_r)
                elif pos.mode=='rec':
                    if pos.side=='l':
                        ret = get_coor_beam_out(pos.aim,pos.i_left,pos.t-pos.tdel,'r',tele_angle=pos.tele_angle_r,beam_angle=pos.beam_angle_r,offset=pos.offset_r)
                    elif pos.side=='r':
                        ret = get_coor_beam_out(pos.aim,pos.i_right,pos.t-pos.tdel,'l',tele_angle=pos.tele_angle_l,beam_angle=pos.beam_angle_l,offset=pos.offset_l)
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
                    if pos.side=='l':
                        ret = get_coor_tele(pos.aim,pos.i_left,pos.t+pos.tdel,'r',tele_angle=pos.tele_angle_r)
                    elif pos.side=='r':
                        ret = get_coor_tele(pos.aim,pos.i_right,pos.t+pos.tdel,'l',tele_angle=pos.tele_angle_l)
                elif pos.mode=='rec':
                    if pos.side=='l':
                        ret = get_coor_tele(pos.aim,pos.i_self,pos.t,'l',tele_angle=pos.tele_angle_l)
                    elif pos.side=='r':
                        ret = get_coor_tele(pos.aim,pos.i_self,pos.t,'r',tele_angle=pos.tele_angle_r)
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
                    if pos.side=='l':
                        ret = get_start_calc(pos.aim,pos.i_self,pos.t+pos.tdel0,'l',pos.tele_angle_l)
                    elif pos.side=='r':
                        ret = get_start_calc(pos.aim,pos.i_self,pos.t+pos.tdel0,'r',pos.tele_angle_r)
                elif pos.mode=='rec':
                    if pos.side=='l':
                        ret = get_start_calc(pos.aim,pos.i_left,pos.t-pos.tdel,'r',pos.tele_angle_r)
                    elif pos.side=='r':
                        ret = get_start_calc(pos.aim,pos.i_right,pos.t-pos.tdel,'l',pos.tele_angle_l)
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
                    if pos.side=='l':
                        ret = get_start_calc(pos.aim,pos.i_left,pos.t+pos.tdel,'r',pos.tele_angle_r)+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
                    elif pos.side=='r':
                        ret = get_start_calc(pos.aim,pos.i_right,pos.t+pos.tdel,'l',pos.tele_angle_l)+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
                elif pos.mode=='rec':
                    if pos.side=='l':
                        ret = get_start_calc(pos.aim,pos.i_self,pos.t-pos.tdel0,'l',pos.tele_angle_l)+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
                    elif pos.side=='r':
                        ret = get_start_calc(pos.aim,pos.i_self,pos.t-pos.tdel0,'r',pos.tele_angle_r)+pos.coor_end[1]*pos.ksi[1]+pos.coor_end[2]*pos.ksi[0]
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos

    def get_direction(self,pos):
        check=False
        while check==False:
            try:
                ret = pos.coor_startbeam[0]
                check=True
            except AttributeError, e:
                self.add_attribute(e,pos)
        return pos








































    def pupil_func(self,f,x,y):
        if (x**2+y**2)**0.5>self.aim.data.D/2.0:
            return np.nan
        else:
            return f(x,y)

    def mean_var(self,i,t,side,ret,mode='mean',Nbins=False,tele_angle_l=False,tele_angle_r=False,beam_angle_l=False,beam_angle_r=False):
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
        #func = lambda x,y: self.pupil_func(func_calc,x,y)

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
        
        func=utils.Object() 
        sampled=utils.Object()

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




        for k in ret:
            if 'adjust' == k:
                for s in side:
                    for i_sel in i:
                        if self.aim.PAAM_deg==1:
                            if s=='l':
                                A = [np.array(getattr(self.aim,'t_adjust')[0][int(i_sel)-1])]
                                A.append(np.array([self.aim.tele_l_ang(i_sel,t) for t in A[0]]))
                            elif s=='r':
                                A = [np.array(getattr(self.aim,'t_adjust')[1][int(i_sel)-1])]
                                A.append(np.array([self.aim.tele_r_ang(i_sel,t) for t in A[0]]))
                        elif self.aim.PAAM_deg==2:
                            A = [np.array(getattr(self.aim,'t_adjust')[s][int(i_sel)])]
                            A.append(np.array(getattr(self.aim,'tele_adjust')[s][int(i_sel)]))
                            
                        B = [A,'value='+k+', mode='+str(mode)]
                        setattr(getattr(getattr(sampled,s),'i'+str(i_sel)),k,B)
            else:
                if option=='both' or option=='function':
                    #setattr(func,k,lambda i,t,side: getattr(self.mean_var(i,t,side,[k],Nbins=Nbins,mode=mode),k))
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
                            #print(getattr(sampled,s),k,B)
                            setattr(getattr(getattr(sampled,s),'i'+str(i_sel)),k,B)
        #self.func = func
        #self.sampled = sampled
        return [func,sampled]


### Calculate values/properties

def get_coor_tele(aim,i,t,side,tele_angle=False):
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

def get_coor_beam_in(aim,i,t,tdel,side,tele_angle_send=False,beam_angle_send=False,tele_angle_rec=False,offset=False,out=3):
    [i_self,i_left,i_right] = utils.i_slr(i)
    check=False
    if aim.data.calc_method=='Abram':
        tdel0 = 0
    elif aim.data.calc_method=='Waluschka':
        tdel0 = tdel
    try:
        if tele_angle_send==False and beam_angle_send==False:
            if side=='l':
                u_not_ab = aim.beam_r_coor(i_left,t-tdel)
            elif side=='r':
                u_not_ab = aim.beam_l_coor(i_right,t-tdel)
            check=True
    except AttributeError:
        check=False
        pass

    if check==False:
        if tele_angle_send is False:
            if side=='l':
                #tele_angle_send = aim.tele_r_ang(i_left,t-tdel)
                tele_angle_send = np.radians(30.0)
            elif side=='r':
                #tele_angle_send = aim.tele_l_ang(i_right,t-tdel)
                tele_angle_send = np.radians(-30.0)
        if beam_angle_send is False:
            if side=='l':
                #beam_angle_send = aim.beam_r_ang(i_left,t-tdel)
                beam_angle_send = 0.0
            elif side=='r':
                #beam_angle_send = aim.beam_l_ang(i_right,t-tdel)
                beam_angle_send = 0.0
        if offset is False:
            if side=='l':
                offset = get_offset(aim,i_left,t-tdel,'r')
            elif side=='r':
                offset = get_offset(aim,i_right,t-tdel,'l')
        elif offset == None:
            offset = 0.0

        if side=='l':
            u_not_ab = methods.beam_coor_out(aim.data,i_left,t-tdel,tele_angle_send,beam_angle_send,offset)[0]
        elif side=='r':
            u_not_ab = methods.beam_coor_out(aim.data,i_right,t-tdel,tele_angle_send,beam_angle_send,offset)[0]

    if aim.data.aberration==False:
        ret = u_not_ab
    elif aim.data.aberration==True:
        if side=='l':
            u_ab = np.linalg.norm(u_not_ab)*(LA.unit(LA.unit(u_not_ab)*aim.data.c+(aim.data.vel.abs(i_self,t-tdel0) - aim.data.vel.abs(i_left,t-tdel))))
        elif side=='r':
            u_ab = np.linalg.norm(u_not_ab)*(LA.unit(LA.unit(u_not_ab)*aim.data.c+(aim.data.vel.abs(i_self,t-tdel0) - aim.data.vel.abs(i_right,t-tdel))))
        if aim.data.relativistic==False:
            ret=u_ab
        else:
            coor = get_coor_tele(aim,i_self,t,side,tele_angle=tele_angle_rec)
            if side=='l':
                velo = (aim.data.vel.abs(i_self,t-tdel0) - aim.data.vel.abs(i_left,t-tdel))
            elif side=='r':
                velo = (aim.data.vel.abs(i_self,t-tdel0) - aim.data.vel.abs(i_right,t-tdel))

            c_vec = LA.unit(u_not_ab)*aim.data.c

            r = coor[0]
            x_prime = LA.unit(velo)
            n_prime = LA.unit(np.cross(velo,r))
            r_prime = LA.unit(np.cross(n_prime,x_prime))

            coor_velo = np.array([r_prime,n_prime,x_prime])
            c_velo = LA.matmul(coor_velo,c_vec)
            v = np.linalg.norm(velo)
            den = 1.0 - ((v/(c**2))*c_velo[2]) #1.0 - ((v/(c**2))*coor_velo[2]) ...adjusted
            num = ((1.0-((v**2)/(c**2)))**0.5)

            ux_prime = (c_velo[2] - v)/den
            ur_prime = (num*c_velo[0])/den
            un_prime = (num*c_velo[1])/den
            c_prime = ux_prime*x_prime + un_prime*n_prime +ur_prime*r_prime
            u_new=LA.unit(c_prime)*np.linalg.norm(u_not_ab)
            
            #ret = np.array([ur_prime,un_prime,ux_prime])
            ret = u_new
#Hier gebleven. Matrix returnen ook bij FAlse relativistic of abberation
    try:
        u_ab
    except NameError:
        u_ab=np.array([np.nan,np.nan,np.nan])

    if out==2:
        ret = [ret,u_not_ab]
    elif out==1:
        ret = [ret]
    elif out==3:
        ret = [ret,u_not_ab,u_ab]


    return ret

def get_coor_beam_out(aim,i,t,side,tele_angle=False,beam_angle=False,offset=False):
    check=False
    
    #print('Test 2: '+str(beam_angle))
    if check==False:
        if tele_angle is False:
            if side == 'l':
                tele_angle = aim.tele_l_ang(i,t)
                #tele_angle = np.radians(-30.0)
            elif side =='r':
                tele_angle = aim.tele_r_ang(i,t)
                #tele_angle = np.radians(30.0)
        elif tele_angle ==None:
            if side == 'l':
                tele_angle = np.radians(-30.0)
            elif side =='r':
                tele_angle = np.radians(30.0)

        if beam_angle is False:
            #print('Test 3: '+str(beam_angle))
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
            #offset = get_offset_tele(aim,i,t,side)
            raise ValueError

        ret = methods.beam_coor_out(aim.data,i,t,tele_angle,beam_angle,offset)

    return ret

def get_offset(aim,i,t,side):
    try:
        ret = aim.offset[side][i](t)
    except TypeError:
        #print(type(side))
        #print(type(i))
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
    try:
        if side=='l':
            ret = aim.tele_l_start(i,t)
        elif side=='r':
            ret = aim.tele_r_start(i,t)
    except AttributeError:
        ret = np.array(aim.data.putp(i,t)) + LA.unit(get_coor_tele(aim,i,t,side,tele_angle=tele_angle)[0])*aim.data.L_tele

    return ret

def values(inp,i,t,side,ksi=[0,0],mode='send',tele_angle_l=False,tele_angle_r=False,beam_angle_l=False,beam_angle_r=False,offset_l=False,offset_r=False,ret=[]):
    [i_self,i_left,i_right] = utils.i_slr(i)
    

    if 'AIM' in str(inp):
        aim = inp
        outp = OUTPUT(aim)
    elif 'OUTPUT' in str(inp):
        aim = inp.aim
        outp = inp
    else:
        raise ValueError
     
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

            #coor_starttele=get_coor_tele(aim,i_self,t,'l',tele_angle=tele_angle_l)
            ##print('Test: '+str(beam_angle_l))
            #coor_startbeam=get_coor_beam_out(aim,i_self,t,'l',tele_angle=tele_angle_l,beam_angle=beam_angle_l,offset=offset_l)
            #[vec_endbeam,vec_endbeam_na,vec_endbeam_ab]=get_coor_beam_in(aim,i_left,t+tdel,tdel,'r',tele_angle_send=tele_angle_l,beam_angle_send=beam_angle_l,tele_angle_rec=tele_angle_r,offset=offset_l)

            #coor_end=get_coor_tele(aim,i_left,t+tdel,'r',tele_angle=tele_angle_r)
            #start=get_start_calc(aim,i_self,t+tdel0,'l',tele_angle_l)
            #end=get_start_calc(aim,i_left,t+tdel,'r',tele_angle_r)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]

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


            #coor_starttele=get_coor_tele(aim,i_self,t,'r',tele_angle=tele_angle_r)
            #coor_startbeam=get_coor_beam_out(aim,i_self,t,'r',tele_angle=tele_angle_r,beam_angle=beam_angle_r,offset=offset_r)       

            #[vec_endbeam,vec_endbeam_na,vec_endbeam_ab]=get_coor_beam_in(aim,i_right,t+tdel,tdel,'l',tele_angle_send=tele_angle_r,beam_angle_send=beam_angle_r,tele_angle_rec=tele_angle_l,offset=offset_r)
            #coor_end=get_coor_tele(aim,i_right,t+tdel,'l',tele_angle=tele_angle_l)

            #start=get_start_calc(aim,i_self,t+tdel0,'r',tele_angle_r)
            #end=get_start_calc(aim,i_right,t+tdel,'l',tele_angle_l)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]

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

            #coor_starttele=get_coor_tele(aim,i_left,t-tdel,'r',tele_angle=tele_angle_r)
            #coor_startbeam=get_coor_beam_out(aim,i_left,t-tdel,'r',tele_angle=tele_angle_r,beam_angle=beam_angle_r,offset=offset_r)
            #[vec_endbeam,vec_endbeam_na,vec_endbeam_ab]=get_coor_beam_in(aim,i_self,t,tdel,'l',tele_angle_send=tele_angle_r,beam_angle_send=beam_angle_r,tele_angle_rec=tele_angle_l,offset=offset_r)
            #coor_end=get_coor_tele(aim,i_self,t,'l',tele_angle=tele_angle_l)
            #
            #start=get_start_calc(aim,i_left,t-tdel,'r',tele_angle_r)
            #end=get_start_calc(aim,i_self,t-tdel0,'l',tele_angle_l)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]
            
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
            
            #coor_starttele=get_coor_tele(aim,i_right,t-tdel,'l',tele_angle=tele_angle_l)
            #coor_startbeam=get_coor_beam_out(aim,i_right,t-tdel,'l',tele_angle=tele_angle_l,beam_angle=beam_angle_l,offset=offset_l)

            #[vec_endbeam,vec_endbeam_na,vec_endbeam_ab]=get_coor_beam_in(aim,i_self,t,tdel,'r',tele_angle_send=tele_angle_l,beam_angle_send=beam_angle_l,tele_angle_rec=tele_angle_r,offset=offset_l)
            #coor_end=get_coor_tele(aim,i_self,t,'r',tele_angle=tele_angle_r)
            #
            #start=get_start_calc(aim,i_right,t-tdel,'l',tele_angle_l)
            #end=get_start_calc(aim,i_self,t-tdel0,'r',tele_angle_r)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]
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
 
    positions=utils.Object()
    positions.method = aim.data.calc_method
    positions.aim=aim
    positions.tele_angle_l = tele_angle_l
    positions.tele_angle_r = tele_angle_r
    positions.beam_angle_l = beam_angle_l
    positions.beam_angle_r = beam_angle_r
    positions.offset_l = offset_l
    positions.offset_r = offset_r
    
    param = ['mode','side','i_self','i_left','i_right','t','ksi','tdel','tdel0','tele_angle_start','tele_angle_end','beam_angle_start','beam_angle_end','aim','offset_start','offset_end']
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

def tele_center_calc(aim,i,t,scale=1,lim=1e-12,max_count=5,print_on=False,value=0,tele_l=False,tele_r=False,beam_l=False,beam_r=False,offset_l=False,offset_r=False):
    [i_self,i_left,i_right] = utils.i_slr(i)
    #print(aim.tele_l_ang,t)
    
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
    
    count=0
    check=False
    while check==False:
        tele_l_old = tele_l
        tele_r_old = tele_r
        pos_send = values(aim,i_self,t,'l',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['off'])
        angx_send = np.sign(pos_send.xoff-value)*abs(np.arctan((pos_send.xoff-value)/pos_send.zoff))
        tele_l = tele_l+angx_send

        pos_rec = values(aim,i_left,t,'r',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['off'])
        angx_rec = np.sign(pos_rec.xoff-value)*abs(np.arctan((pos_rec.xoff-value)/pos_rec.zoff))
        #print('Check 1: '+ str(pos_send.xoff)+ ', '+str(pos_rec.xoff))
        tele_r = tele_r+angx_rec

        count=count+1

        if count>=max_count:
            check=True
            mode='Maximum counts'
        
        if tele_l_old-lim<=tele_l and tele_l<=tele_l_old+lim and tele_r_old-lim<=tele_r and tele_r<=tele_r_old+lim:
            chek=True
            mode='Converged'
        if abs(angx_send)<=lim and abs(angx_rec)<=lim:
            chek=True
            mode='Converged'
        if print_on:
            print(tele_l,tele_r)
    #print(pos_rec.xoff)

    return [[tele_l,tele_r], mode]

def PAAM_center_calc(aim,i,t,para='yoff',scale=1,lim=1e-12,max_count=5,print_on=False,tele_l=None,tele_r=None,beam_l=None,beam_r=None,offset_l=None,offset_r=None,value=0,method='iter',margin=0.01):
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
    
    if method=='iter': 
        count=0
        check=False
        while check==False:
            beam_l_old = beam_l
            beam_r_old = beam_r
            pos_send = values(aim,i_self,t,'l',tele_angle_l=tele_l,tele_angle_r=tele_l,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['off'])
            angy_send = np.sign(pos_send.yoff)*abs(np.arctan(pos_send.yoff/pos_send.zoff))
            beam_l = beam_l-angy_send

            pos_rec = values(aim,i_left,t,'r',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['off'])
            angy_rec = np.sign(pos_rec.yoff)*abs(np.arctan(pos_rec.yoff/pos_rec.zoff))
            beam_r = beam_r-angy_rec

            count=count+1

            if count>=max_count:
                check=True
                mode='Maximum counts'
            if beam_l_old-lim<=beam_l and beam_l<=beam_l_old+lim and beam_r_old-lim<=beam_r and beam_r<=beam_r_old+lim:
                chek=True
                mode='Converged'
            if abs(angy_send)<=lim and abs(angy_rec)<=lim:
                chek=True
                mode='Converged'
            if print_on:
                print(beam_l,beam_r)
        #print(pos_rec.yoff)

    elif method=='solve':
        solve_l = lambda beam_l_solve: getattr(values(aim,i_self,t,'l',tele_angle_l=tele_l,tele_angle_r=tele_l,beam_angle_l=beam_l_solve,beam_angle_r=beam_r,ret=[para]),para) - value
        beam_l = scipy.optimize.brentq(solve_l,-margin,margin,xtol=lim)
        solve_r = lambda beam_r_solve: getattr(values(aim,i_left,t,'r',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r_solve,ret=[para]),para) - value
        beam_r = scipy.optimize.brentq(solve_r,-margin,margin,xtol=lim)
        mode='Converged'

    return [[beam_l,beam_r], mode]



def PAAM_wavefront_calc(aim,i,t,side,lim=1e-12):
    if side=='l':
        angy = lambda beam: values(aim,i,t,'l',mode='send',beam_angle_l=beam,ret=['angy_wf_send']).angy_wf_send
    elif side=='r':
        angy = lambda beam: values(aim,i,t,'r',mode='send',beam_angle_r=beam,ret=['angy_wf_send']).angy_wf_send

    try:
        ret = scipy.optimize.brentq(angy,-1e-5,1e-5,xtol=lim)
    except ValueError,e:
        if str(e)=='f(a) and f(b) must have different signs':
            ret=np.nan
    return ret

def tele_wavefront_calc(aim,i_l,t,method,scale=1,lim=1e-12,max_count=20,tele_angle_l=None,tele_angle_r=None,beam_l=None,beam_r=None,offset_l=False,offset_r=False,print_on=False): #The sending is always the left telescope and the receiving the right one
    i_r = utils.i_slr(i_l)[1]
    
    tdel = aim.data.L_sl_func_tot(i_l,t)
    t_l=t
    t_r=t+tdel #...because considering traveling time makes it more complex (solve for functions)
    #t_l=t
    #t_r=t #...because considering traveling time makes it more complex (solve for functions)

    if tele_angle_l is False:
        tele_angle_l=aim.tele_l_ang(i_l,t_l)
    elif tele_angle_l == None:
        tele_angle_l=np.radians(-30.0)
    if tele_angle_r is False:
        tele_angle_r=aim.tele_r_ang(i_r,t_r)
    elif tele_angle_r ==None:
        tele_angle_r = np.radians(30.0)

    calc_l=100
    calc_r=100
    calc_l_old=calc_l
    calc_r_old=calc_r

    count=0
    lim_val=100

    para = 'angx_wf_send'
    
    tdel = aim.data.L_sl_func_tot(i_l,t)

    if method=='iter':
        while count<max_count:
            calc_l_old = calc_l
            calc_r_old = calc_r
            tele_angle_l_old = tele_angle_l
            tele_angle_r_old = tele_angle_r
            pos_l = values(aim,i_l,t_l,'l',mode='send',tele_angle_l=tele_angle_l,tele_angle_r=tele_angle_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r)
            calc_l = getattr(getattr(OUTPUT(aim),'get_'+para)(pos_l),para)
            tele_angle_r = tele_angle_r-scale*calc_l

            pos_r = values(aim,i_r,t_r+tdel,'r',mode='send',tele_angle_l=tele_angle_l,tele_angle_r=tele_angle_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r)
            calc_r = getattr(getattr(OUTPUT(aim),'get_'+para)(pos_r),para)
            tele_angle_l = tele_angle_l-scale*calc_r

            count=count+1
            lim_val = max(abs(abs(calc_l_old)-abs(calc_l)),abs(abs(calc_r_old)-abs(calc_r)))
            if print_on:
                print(tele_angle_l,tele_angle_r)
                print(calc_l,calc_r)
                print(max(abs(calc_l),abs(calc_r)),max(abs(calc_l),abs(calc_r))<lim)
                print(lim_val)
                print('')

            #if abs(calc_l_old)-lim<=abs(calc_l) and abs(calc_l)<=abs(calc_l_old)+lim and abs(calc_r_old)-lim<=abs(calc_r) and abs(calc_r)<=abs(calc_r_old)+lim:
            if lim_val<=lim:
                mode = 'Result is converged'
                if print_on:
                    print(mode)
                break

            if count>= max_count:
                mode = 'Maximum iteration limit has been reached'
                if print_on:
                    print(mode)
                break
    return [[tele_angle_l,tele_angle_r],mode]

def get_tele_wavefront(aim,i,t,side,method,scale=1,lim=1e-12,max_count=20,print_on=False,value=0.0,tele_angle_l=None,tele_angle_r=None,beam_l=None,beam_r=None,offset_l=False,offset_r=False): 
    if side=='l':
        i_l = i
        tdel=0
    elif side=='r':
        i_l = utils.i_slr(i)[2]
        i_r = i
        tdel = aim.data.L_rr_func_tot(i_r,t)
    
    try:
        ang = tele_wavefront_calc(aim.aim0,i_l,t-tdel,aim.aimset.tele_method_solve,lim=aim.aimset.limit_angx,scale=scale,max_count=max_count,print_on=print_on,tele_angle_l=tele_angle_l,tele_angle_r=tele_angle_r,beam_l=beam_l,beam_r=beam_r,offset_l=offset_l,offset_r=offset_r)
        aim_sel = aim.aim0
    except AttributeError:
        ang = tele_wavefront_calc(aim,i_l,t-tdel,aim.aimset.tele_method_solve,lim=aim.aimset.limit_angx,scale=scale,max_count=max_count,print_on=print_on,tele_angle_l=tele_angle_l,tele_angle_r=tele_angle_r,beam_l=beam_l,beam_r=beam_r,offset_l=offset_l,offset_r=offset_r)
        aim_sel=aim


    if side=='l':
        pos = values(aim_sel,i_l,t,'l',ksi=[0,0],mode='send',tele_angle_l=ang[0][0],tele_angle_r=ang[0][1],beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['angx_wf_send']).angx_wf_send
        ret=ang[0][0]
    elif side=='r':
        pos = values(aim_sel,i_r,t,'r',ksi=[0,0],mode='send',tele_angle_l=ang[0][0],tele_angle_r=ang[0][1],beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['angx_wf_send']).angx_wf_send
        ret=ang[0][1]
    
    if print_on:
        print(ret,pos)

    return ret









#def tele_wavefront_calc_old(aim,i_send,t,para,method,scale=1,lim=1e-12,max_count=5,print_on=False): #The sending is always the left telescope and the receiving the right one
#    # Optimization per 8 seconds accurate
#    
#    i_rec = utils.i_slr(i_send)[1]
#
#    if i_send==0:
#        i_send=3
#    if i_rec==0:
#        i_rec=3
#    t_send=t
#    t_rec=t_send #...because considering traveling time makes it more complex (solve for functions)
#
#    tele_angle_l=aim.tele_l_ang(i_send,t)
#    tele_angle_r=aim.tele_r_ang(i_rec,t)
#    
#    calc_send=100
#    calc_rec=100
#    calc_send_old=calc_send
#    calc_rec_old=calc_send
#    
#    count=0
#    lim_val=100
#    
#    if method=='iter':
#        while count<max_count:
#            calc_send_old = calc_send
#            calc_rec_old = calc_rec
#            tele_angle_l_old = tele_angle_l
#            tele_angle_r_old = tele_angle_r
#
#            pos_send = values(aim,i_send,t_send,'l',mode='send',tele_angle_l=tele_angle_l,tele_angle_r=tele_angle_r)
#            calc_rec = getattr(getattr(OUTPUT(aim),'get_'+para)(pos_send),para)
#            tele_angle_r = tele_angle_r-scale*calc_rec
#            pos_rec = values(aim,i_rec,t_rec,'r',mode='send',tele_angle_l=tele_angle_l,tele_angle_r=tele_angle_r)
#            calc_send = getattr(getattr(OUTPUT(aim),'get_'+para)(pos_rec),para)
#            tele_angle_l = tele_angle_l-scale*calc_send
#            count=count+1
#            
#            lim_val = max(abs(calc_send),abs(calc_rec))
#            
#            if print_on:
#                print(tele_angle_l,tele_angle_r)
#                print(calc_send,calc_rec)
#                print(max(abs(calc_send),abs(calc_rec)),max(abs(calc_send),abs(calc_rec))>lim)
#                print('')
#            
#            if abs(calc_send_old)-lim<=abs(calc_send) and abs(calc_send)<=abs(calc_send_old)+lim and abs(calc_rec_old)-lim<=abs(calc_rec) and abs(calc_rec)<=abs(calc_rec_old)+lim:
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
#
#    elif method=='solve':
#        mode=method
#        tele_angle_l=aim.tele_l_ang(i_send,t)
#        tele_angle_r=aim.tele_r_ang(i_rec,t)
#        lim_val=100
#
#        count=0
#        val_send_old = tele_angle_l
#        val_rec_old = tele_angle_r
#        pos_send = lambda tele_l: values(aim,i_send,t_send,'l',mode='send',tele_angle_l=tele_l,tele_angle_r=tele_angle_r)
#        calc_send = lambda tele_l: getattr(getattr(OUTPUT(aim),'get_'+para)(pos_send(tele_l)),para)
#        tele_angle_l = scipy.optimize.brentq(calc_send,np.radians(-30-5),np.radians(-30+5),xtol=lim)
#        pos_rec = lambda tele_r: values(aim,i_rec,t_rec,'r',mode='send',tele_angle_l=tele_angle_l,tele_angle_r=tele_r)
#        calc_rec =lambda tele_r:  getattr(getattr(OUTPUT(aim),'get_'+para)(pos_rec(tele_r)),para)
#        tele_angle_r = scipy.optimize.brentq(calc_rec,np.radians(30-5),np.radians(30+5),xtol=lim)
#         
#        lim_val=max(abs(calc_send(tele_angle_l)),abs(calc_rec(tele_angle_r)))
#        
#        if print_on:
#            print(tele_angle_l,tele_angle_r)
#            print(calc_send(tele_angle_l),calc_rec(tele_angle_r))
#            print('')
#
#        #count=count+1
#        #if count>=max_count:
#        #    mode='Maximum iterations reached'
#        #    if print_on:
#        #        print(mode)
#        #    break
#
#        #if val_send_old==tele_angle_l and val_rec_old==tele_angle_r:
#        #    mode='Convergence received'
#        #    if print_on:
#        #        print(mode)
#        #    break
#
#    return [[tele_angle_l,tele_angle_r],mode]
