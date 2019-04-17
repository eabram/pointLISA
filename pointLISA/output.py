from imports import *
import inspect

class OUTPUT():
    def __init__(self,aim,**kwargs):
        self.aim = aim

        ### Obtain parameters
        for k, value in parameters.__dict__.items():
            if '__' not in k:
                globals()[k] = value
                setattr(self,k,value)
    
    def pupil(self,Nbins=2): #Aperture
        D_calc=self.D

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



    def add_attribute(self,e,pos):
        func = 'get_'+str(e).split("attribute '")[-1].replace("'","")
        pos = getattr(self,func)(pos)

        return pos
    
    def get_xoff(self,pos):
        setattr(pos,inspect.stack()[0][3].split('get_')[1],LA.matmul(pos.coor_start,pos.end-pos.start)[2])
        return pos

    def get_yoff(self,pos):
        setattr(pos,inspect.stack()[0][3].split('get_')[1],LA.matmul(pos.coor_start,pos.end-pos.start)[1])
        return pos
    
    def get_zoff(self,pos):
        setattr(pos,inspect.stack()[0][3].split('get_')[1],LA.matmul(pos.coor_start,pos.end-pos.start)[0])
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


    
    def get_R(self,pos,precision=0):
        check=False
        while check==False:
            try:
                if precision==0:
                    ret=pos.zoff
                elif precision==1:
                    try:
                       [piston,z_extra] = wfe.z_solve(xoff,yoff,zoff,ret='all')
                    except:
                        [piston,z_extra] = [np.nan,np.nan]
                    ret = self.R(piston)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos) 
        return pos

    def get_R_vec(self,pos):
        check=False
        while check==False:
            try:
                ret = np.array([(pos.R**2-pos.xoff**2-pos.yoff**2)**0.5,pos.yoff,pos.xoff])
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
                ret = LA.matmul(np.linalg.inv(pos.coor_start),pos.R_vec)
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
                ret = LA.matmul(pos.coor_end,-pos.R_vec_origin)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos


  
    def get_tele_vec(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_start,-pos.coor_end[0])
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
                ret = np.sign(pos.R_vec[2])*abs(np.arctan(pos.R_vec[2]/pos.R_vec[0]))
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
                ret = np.sign(pos.R_vec[1])*abs(np.arctan(pos.R_vec[1]/pos.R_vec[0]))
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

    def get_angx_func_rec(self,pos):
        check=False
        while check==False:
            try:
                #ret = (pos.angx_tele-pos.angx_R)
                ret = np.arctan(abs(pos.R_vec_tele_rec[2]/pos.R_vec_tele_rec[0]))*np.sign(pos.R_vec_tele_rec[2])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_angy_func_rec(self,pos):
        check=False
        while check==False:
            try:
                #ret = (pos.angy_tele-pos.angy_R)
                ret = np.arctan(abs(pos.R_vec_tele_rec[1]/pos.R_vec_tele_rec[0]))*np.sign(pos.R_vec_tele_rec[1])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos
    
    def get_bd_original_frame(self,pos):
        check=False
        while check==False:
            try:
                ret = np.array(pos.coor_start[0])
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_bd_receiving_frame(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.matmul(pos.coor_end,pos.bd_original_frame)
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_piston(self,pos):
        check=False
        while check==False:
            try:
                ret = self.z_solve(pos.xoff,pos.yoff,pos.zoff,ret='all')[0]
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
                ret = self.z_solve(pos.xoff,pos.yoff,pos.zoff,ret='all')[1]
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
                ret = np.arccos(-pos.bd_receiving_frame[0]/np.linalg.norm(pos.bd_receiving_frame))
                setattr(pos,inspect.stack()[0][3].split('get_')[1],ret)
                check=True
            except AttributeError,e:
                #print(e)
                self.add_attribute(e,pos)
        return pos

    def get_FOV_wavefront(self,pos):
        check=False
        while check==False:
            try:
                ret = LA.angle(-pos.R_vec_origin,pos.coor_end[0])
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












    def values(self,i,t,side,ksi=[0,0],mode='send',angles=False,ret=[]):
        [i_self,i_left,i_right] = utils.i_slr(i)

        if mode=='send':
            if side=='l':
                tdel = self.aim.data.L_sl_func_tot(i_self,t)
                if self.aim.data.calc_method=='Waluschka':
                    tdel0=tdel
                elif self.aim.data.calc_method=='Abram':
                    tdel0=0
                tele_ang=self.aim.tele_l_ang(i_self,t+tdel0)
                coor_start=self.aim.tele_l_coor(i_self,t)
                coor_end=self.aim.tele_r_coor(i_left,t+tdel)
                start=self.aim.tele_l_start(i_self,t+tdel0)
                end=self.aim.tele_r_start(i_left,t+tdel)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]

            elif side=='r':
                tdel = self.aim.data.L_sr_func_tot(i_self,t)
                if self.aim.data.calc_method=='Waluschka':
                    tdel0=tdel
                elif self.aim.data.calc_method=='Abram':
                    tdel0=0
                tele_ang=self.aim.tele_r_ang(i_self,t+tdel0)
                coor_start=self.aim.tele_r_coor(i_self,t)
                coor_end=self.aim.tele_l_coor(i_right,t+tdel)
                start=self.aim.tele_r_start(i_self,t+tdel0)
                end=self.aim.tele_l_start(i_right,t+tdel)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]

        elif mode=='rec':
            if side=='l':
                tdel = self.aim.data.L_rl_func_tot(i_self,t)
                if self.aim.data.calc_method=='Waluschka':
                    tdel0=tdel
                elif self.aim.data.calc_method=='Abram':
                    tdel0=0
                tele_ang=self.aim.tele_r_ang(i_left,t-tdel)
                coor_start=self.aim.tele_r_coor(i_left,t-tdel)
                coor_end=self.aim.tele_l_coor(i_self,t)
                start=self.aim.tele_r_start(i_left,t-tdel)
                end=self.aim.tele_l_start(i_self,t-tdel0)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]

            elif side=='r':
                tdel = self.aim.data.L_rr_func_tot(i_self,t)
                if self.aim.data.calc_method=='Waluschka':
                    tdel0=tdel
                elif self.aim.data.calc_method=='Abram':
                    tdel0=0
                tele_ang=self.aim.tele_l_ang(i_right,t-tdel)
                coor_start=self.aim.tele_l_coor(i_right,t-tdel)
                coor_end=self.aim.tele_r_coor(i_self,t)
                start=self.aim.tele_l_start(i_right,t-tdel)
                end=self.aim.tele_r_start(i_self,t-tdel0)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]
        
        positions=utils.Object()
        positions.start=start
        positions.end=end
        positions.coor_start=coor_start
        positions.coor_end=coor_end
        positions.tdel=tdel
        positions.tdel0=tdel0
        positions.ksi=ksi

        for r in ret:
            #print(r)
            if r not in positions.__dict__.keys():
                positions_new = getattr(self,'get_'+r)(positions)
                del positions
                positions = positions_new

        return positions
        

                









