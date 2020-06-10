from pointLISA import *

# in the CONSTELLATION class an ORBIT object is used which containes the imported or constructed orbital coordinates of the spacecrafts. In CONSTELLATION different vectors are being constructed and calculaded (r, n, v, u, L) as well as angles (PAA, breathing)

class CONSTELLATION():
    def __init__(self,settings=None,orbit_file=None,**kwargs):
        print('##############################')
        input_param = const.get_settings(settings_input=settings,select='constellationset',kwargs=kwargs)
        setattr(input_param,'filename',orbit_file)
        self.constellationset = input_param
        self.param = self.get_parameters()
        
        self.PAA_func()
        self.settings = settings
        print('##############################')
        print('')

    def get_parameters(self):
        param = utils.Object()
        for k in parameters_all.keys():
            setattr(param,k,parameters_all[k])

        return param

    def get_PAA(self,i,t,side,give='inp'):
        ''' Obtains the Point-ahead angles for the components inp (inplane), out (outplane) and tot (total)'''
        if side=='l':
            v1 = self.v_l(i,t)
            v2 = -self.u_l(i,t)
        elif side=='r':
            v1 = self.v_r(i,t)
            v2 = -self.u_r(i,t)
        n = self.n_func(i,t)
        r = self.r_func(i,t)

        return LA.ang_in_out_tot(v1,v2,n,r,give=give)

    # Velocity
    def velocity_abs_calc(self,i_select,t,hstep):
        '''Returns the velocity vector of a spacecraft'''
        v = (np.array(self.putp(i_select,np.float64(t+hstep)))-np.array(self.putp(i_select,t)))/hstep
        return v

    def velocity_abs(self):
        '''Returns the velocity vector in a function'''
        hstep = np.float64(self.constellationset.hstep)
        v_ret = lambda i,time: self.velocity_abs_calc(i,time,hstep)
        try:
            self.vel
        except AttributeError:
            self.vel = utils.Object()
        self.vel.abs = v_ret

        return self.vel.abs

    def velocity_calc(self,i,time,hstep,side,rs):
        '''Calculates the velocity components'''
        [i_OBJ,i_next] = utils.const.i_slr(i,side=side)
        v_pos_l = self.v_l_stat_func_tot(i_OBJ,time)
        v_pos_r = self.v_r_stat_func_tot(i_OBJ,time)
        n = np.cross(v_pos_r,v_pos_l)
        n = LA.unit(n)
        if side=='l':
            v_pos = v_pos_l
        elif side=='r':
            v_pos = v_pos_r

        pos_OBJ = np.array(self.putp(i_OBJ,time))
        pos_next = np.array(self.putp(i_next,time))
        pos_OBJ_h = np.array(self.putp(i_OBJ,time+np.float64(hstep)))
        pos_next_h = np.array(self.putp(i_next,time+np.float64(hstep)))
        v = ((pos_next_h-pos_next) - (pos_OBJ_h - pos_OBJ))/hstep

        v_out = n*(np.dot(v,n))
        v_arm = v*(np.dot(LA.unit(v),LA.unit(v_pos)))
        v_in = v - v_out - v_arm

        v_out_sign = np.sign(np.dot(v_out,n))
        v_arm_sign = np.sign(np.dot(LA.unit(v),LA.unit(v_pos)))
        v_in_sign = np.sign(np.dot(v_in,v_pos_r - v_pos_l))

        v_out_mag = np.linalg.norm(v_out)*v_out_sign
        v_in_mag = np.linalg.norm(v_in)*v_in_sign
        v_arm_mag = np.linalg.norm(v_arm)*v_arm_sign

        ret =  [v,v_in,v_out,v_arm,v_in_mag,v_out_mag,v_arm_mag]

        return ret[rs]

    def velocity_func(self,hstep=1.0):
        '''Returns functions of al velocity components'''
        hstep = np.float64(hstep)

        try:
            self.vel
        except AttributeError:
            self.vel = Object()
        self.vel.l= lambda i,time: self.velocity_calc(i,time,hstep,'l',0)
        self.vel.in_l = lambda i,time: self.velocity_calc(i,time,hstep,'l',1)
        self.vel.out_l = lambda i,time: self.velocity_calc(i,time,hstep,'l',2)
        self.vel.arm_l = lambda i,time: self.velocity_calc(i,time,hstep,'l',3)
        self.vel.in_mag_l = lambda i,time: self.velocity_calc(i,time,hstep,'l',4)
        self.vel.out_mag_l = lambda i,time: self.velocity_calc(i,time,hstep,'l',5)
        self.vel.arm_mag_l = lambda i,time: self.velocity_calc(i,time,hstep,'l',6)
        self.vel.r= lambda i,time: self.velocity_calc(i,time,hstep,'r',0)
        self.vel.in_r = lambda i,time: self.velocity_calc(i,time,hstep,'r',1)
        self.vel.out_r = lambda i,time: self.velocity_calc(i,time,hstep,'r',2)
        self.vel.arm_r = lambda i,time: self.velocity_calc(i,time,hstep,'r',3)
        self.vel.in_mag_r = lambda i,time: self.velocity_calc(i,time,hstep,'r',4)
        self.vel.out_mag_r = lambda i,time: self.velocity_calc(i,time,hstep,'r',5)
        self.vel.arm_mag_r = lambda i,time: self.velocity_calc(i,time,hstep,'r',6)

        return 0

    def get_armvec_func(self,i,t,side):
        '''Obtains the functions of the distance vectors between two spacecrafts'''
        [i_self,i_left,i_right] = const.i_slr(i)
        if side=='l':
            i_next = i_left
        elif side=='r':
            i_next = i_right
        
        return self.putp(i_next,t) - self.putp(i_self,t)

    def func_over_sc(self,func_tot):
        '''Makes from a list of funcions a function (wih two variables)'''
        f = lambda i,t: func_tot[i-1](t)

        return f

    def PAA_func(self):
        '''Obtains functions of vectors and angles (PAA, brething, n, r, u, v, L'''
        print('Importing Orbit')
        tic=time.clock()
        Orbit=ORBIT.ORBIT(input_param=self.constellationset.__dict__)
        try:
            print(str(Orbit.linecount)+' datapoints')
        except AttributeError:
            print('Obtained orbital function')
        self.orbit = Orbit
        self.LISA = self.orbit.LISA
        self.putp = self.orbit.putp
        print('')
        print('Done in '+str(time.clock()-tic)+' seconds')
        self.SC = range(1,4) 

        # Calculations       
        self.L_sl = lambda i,t: self.send_func_new(i,'l',mode='send',give='L')(t)
        self.L_sr = lambda i,t: self.send_func_new(i,'r',mode='send',give='L')(t)
        self.L_rl = lambda i,t: self.send_func_new(i,'l',mode='rec',give='L')(t)
        self.L_rr = lambda i,t: self.send_func_new(i,'r',mode='rec',give='L')(t)
        self.v_l = lambda i,t: self.send_func_new(i,'l',mode='send',give='v')(t)
        self.v_r = lambda i,t: self.send_func_new(i,'r',mode='send',give='v')(t)
        self.u_l = lambda i,t: self.send_func_new(i,'l',mode='rec',give='v')(t)
        self.u_r = lambda i,t: self.send_func_new(i,'r',mode='rec',give='v')(t)
 
        self.v_l_stat = lambda i,t: self.get_armvec_func(i,t,'l')
        self.v_r_stat = lambda i,t: self.get_armvec_func(i,t,'r')

        self.n_func = lambda i,t: LA.unit(np.cross(self.v_l_stat(i,t),self.v_r_stat(i,t)))
        self.r_func = lambda i,t: (self.v_l_stat(i,t)+self.v_r_stat(i,t))/2.0

        self.v_l_in = lambda i,t: LA.inplane_outplane(self.v_l(i,t),self.n_func(i,t))[0]
        self.v_r_in = lambda i,t: LA.inplane_outplane(self.v_r(i,t),self.n_func(i,t))[0]
        self.u_l_in = lambda i,t: LA.inplane_outplane(self.u_l(i,t),self.n_func(i,t))[0]
        self.u_r_in = lambda i,t: LA.inplane_outplane(self.u_r(i,t),self.n_func(i,t))[0]
        self.v_l_out = lambda i,t: LA.inplane_outplane(self.v_l(i,t),self.n_func(i,t))[1]
        self.v_r_out = lambda i,t: LA.inplane_outplane(self.v_r(i,t),self.n_func(i,t))[1]
        self.u_l_out = lambda i,t: LA.inplane_outplane(self.u_l(i,t),self.n_func(i,t))[1]
        self.u_r_out = lambda i,t: LA.inplane_outplane(self.u_r(i,t),self.n_func(i,t))[1]

        #--- Obtaining PAA --- 
        PAA = utils.Object()
        PAA.inp =  lambda i,t,s: self.get_PAA(i,t,s,give='inp')
        PAA.out =  lambda i,t,s: self.get_PAA(i,t,s,give='out')
        PAA.tot =  lambda i,t,s: self.get_PAA(i,t,s,give='tot')
        self.PAA = PAA
        
        #--- Obtaining breathing angles ---
        self.ang_breathing_din = lambda i, time: LA.angle(self.v_l(i,time),self.v_r(i,time))
        self.ang_breathing_in = lambda i, time: LA.angle(self.u_l(i,time),self.u_r(i,time))
        self.ang_breathing_stat = lambda i, time: LA.angle(self.v_l_stat(i,time),self.v_r_stat(i,time))
        
        #--- Obtaining Velocity
        self.velocity_abs()
        self.velocity_func()

        try:
            self.t_all
        except AttributeError:
            self.t_all = self.orbit.t
        if 'int' in str(type(self.constellationset.length_calc)):
            if self.t_all[-1]>(self.constellationset.length_calc+1)*self.param.day2sec:
                loc = calc.get_nearest_smaller_value(self.orbit.t,self.constellationset.length_calc*self.param.day2sec)
                self.t_all = self.orbit.t[0:loc+1]


        return self

    def solve_L_PAA(self,t,pos_self,pos_left,pos_right,select='sl',i=False):
        '''Calculate the photon traveling time along one of the six laserlinks'''
        t_guess = np.linalg.norm(np.array(self.putp(1,0)) - np.array(self.putp(2,0)))/c

        if select=='sl' or select=='rl':
            s1 = lambda x: pos_left(x)
        elif select=='sr' or select=='rr':
            s1 = lambda x: pos_right(x)

        s2 = lambda x: pos_self(x)
        x_0 = t
        if select=='sl' or select=='sr':
            s3 = lambda dt: s1(x_0+dt) - s2(x_0)
        elif select=='rl' or select=='rr':
            s3 = lambda dt: -s1(x_0-dt) + s2(x_0)
        s4 = lambda dt: np.linalg.norm(s3(dt))
        s5 = lambda dt: s4(dt) - c*dt

        try:
            res = scipy.optimize.brentq(s5,0,t_guess*4)
        except ValueError,e:
            if str(e)=='f(a) and f(b) must have different signs':
                res=np.nan

        return res

    def L_PAA(self,pos_self,pos_left,pos_right,i=False):
        '''Obtain time of flight of beam between spacecrafts'''
        selections = ['sl','sr','rl','rr']

        L_sl_func =  lambda time: self.solve_L_PAA(time,pos_self,pos_left,pos_right,select=selections[0],i=i)
        L_sr_func =  lambda time: self.solve_L_PAA(time,pos_self,pos_left,pos_right,select=selections[1],i=i)
        L_rl_func =  lambda time: self.solve_L_PAA(time,pos_self,pos_left,pos_right,select=selections[2],i=i)
        L_rr_func =  lambda time: self.solve_L_PAA(time,pos_self,pos_left,pos_right,select=selections[3],i=i)

        return [L_sl_func,L_sr_func,L_rl_func,L_rr_func]

    def relativistic_aberrations(self,i,t,u):
        '''Adjust vecor u by adding the angle caused by aberration'''
        relativistic=self.constellationset.relativistic

        if relativistic==True:
            V = -self.vel.abs(i,t)
            V_mag = np.linalg.norm(V)
            u_mag = np.linalg.norm(u)
            c_vec = LA.unit(u)*c

            velo = V
            coor = calc.coor_SC(self,i,t)
            r=coor[0]
            x_prime = LA.unit(velo)
            n_prime = LA.unit(np.cross(velo,r))
            r_prime = LA.unit(np.cross(n_prime,x_prime))
            coor_velo = np.array([r_prime,n_prime,x_prime])
            c_velo = LA.matmul(coor_velo,c_vec)
            v = np.linalg.norm(velo)
            den = 1.0 - ((v/(c**2))*coor_velo[2])
            num = ((1.0-((v**2)/(c**2)))**0.5)

            ux_prime = (c_velo[2] + v)/den
            ur_prime = (num*c_velo[0])/den
            un_prime = (num*c_velo[1])/den
            c_prime = ux_prime*x_prime + un_prime*n_prime +ur_prime*r_prime
            u_new = LA.unit(c_prime)*u_mag

        elif relativistic==False:
            V = -self.vel.abs(i,t)
            u_mag = np.linalg.norm(u)
            c_vec = LA.unit(u)*c

            u_new = LA.unit(c_vec+V)*u_mag
        else:
            print('Error')

            return u_new

    def send_func(self,i):
        '''Uses previous defined functions to return the vecors L, u, v, r and n'''
        [i_self,i_left,i_right] = utils.const.i_slr(i)

        pos_left = lambda time: self.putp(i_left,time)
        pos_self = lambda time: self.putp(i_self,time)
        pos_right = lambda time: self.putp(i_right,time)

        [L_sl,L_sr,L_rl,L_rr] = self.L_PAA(pos_self,pos_left,pos_right,i=i_self)

        #Abram2018
        v_send_l0 = lambda t: pos_left(t+L_sl(t)) - pos_self(t)
        v_send_r0 = lambda t: pos_right(t+L_sr(t)) - pos_self(t)
        v_rec_l0 = lambda t: pos_self(t) - pos_left(t - L_rl(t))
        v_rec_r0 = lambda t: pos_self(t) - pos_right(t - L_rr(t))

        if self.constellationset.aberration==False:
            v_send_l = v_send_l0
            v_send_r = v_send_r0
            v_rec_l = v_rec_l0
            v_rec_r = v_rec_r0
        elif self.constellationset.aberration==True:
            v_send_l = lambda t: self.relativistic_aberrations(i,t,v_send_l0(t))
            v_send_r = lambda t: self.relativistic_aberrations(i,t,v_send_r0(t))
            v_rec_l = lambda t: self.relativistic_aberrations(i,t,v_rec_l0(t))
            v_rec_r = lambda t: self.relativistic_aberrations(i,t,v_rec_r0(t))

        return [[v_send_l,v_send_r,v_rec_l,v_rec_r],[L_sl,L_sr,L_rl,L_rr]]

    def send_func_new(self,i,side,mode='send',give='L'):
        '''Uses previous defined functions to return the vecors L, u, v, r and n'''
        [i_self,i_left,i_right] = utils.const.i_slr(i)

        pos_left = lambda time: self.putp(i_left,time)
        pos_self = lambda time: self.putp(i_self,time)
        pos_right = lambda time: self.putp(i_right,time)

        [L_sl,L_sr,L_rl,L_rr] = self.L_PAA(pos_self,pos_left,pos_right,i=i_self)
        
        if mode=='send':
            if side=='l':
                if give=='L':
                    ret = L_sl
                elif give=='v':
                    v_send_l0 = lambda t: pos_left(t+L_sl(t)) - pos_self(t)
                    ret = lambda t: self.relativistic_aberrations(i,t,v_send_l0(t))

            elif side=='r':
                if give=='L':
                    ret = L_sr
                elif give=='v':
                    v_send_r0 = lambda t: pos_right(t+L_sr(t)) - pos_self(t)

                    ret = lambda t: self.relativistic_aberrations(i,t,v_send_r0(t))

        if mode=='rec':
            if side=='l':
                if give=='L':
                    ret = L_rl
                elif give=='v':
                    v_rec_l0 = lambda t: pos_self(t) - pos_left(t - L_rl(t))
                    ret = v_rec_l = lambda t: self.relativistic_aberrations(i,t,v_rec_l0(t))
            elif side=='r':
                if give=='L':
                    ret = L_rr
                elif give=='v':
                    v_rec_r0 = lambda t: pos_self(t) - pos_right(t - L_rr(t))
                    ret = lambda t: self.relativistic_aberrations(i,t,v_rec_r0(t))

        return ret
