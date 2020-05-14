from pointLISA import *

# in the STAT class an ORBIT object is used which containes the imported or constructed orbital coordinates of the spacecrafts. In STAT different vectors are being constructed and calculaded (r, n, v, u, L) as well as angles (PAA, breathing)

class STAT():
    def __init__(self,settings=None,orbit_file=None,**kwargs):
        input_param = const.get_settings(settings_input=settings,select='stat',kwargs=kwargs)
        setattr(input_param,'filename',orbit_file)
        self.stat = input_param
        self.param = self.get_parameters()
        
        self.PAA_func()
        self.setings = settings

    def get_parameters(self):
        param = utils.Object()
        for k in parameters_all.keys():
            setattr(param,k,parameters_all[k])

        return param

















    #PAA angles
    def calc_PAA_ltot(self,i,t):
        '''Returns the total PAA for the left telecope'''
        calc_ang=LA.angle(self.v_l_func_tot(i,t),-self.u_l_func_tot(i,t))
        return calc_ang

    def calc_PAA_lin(self,i,t):
        '''Returns the inplane PAA for the left telecope'''
        calc_ang=LA.ang_in_out(self.v_l_func_tot(i,t),-self.u_l_func_tot(i,t),self.n_func(i,t),self.r_func(i,t),give='in')
        return calc_ang

    def calc_PAA_lout(self,i,t):
        '''Returns the outplanr PAA for the left telecope'''
        calc_ang=LA.ang_in_out(self.v_l_func_tot(i,t),-self.u_l_func_tot(i,t),self.n_func(i,t),self.r_func(i,t),give='out')
        return calc_ang

    def calc_PAA_rtot(self,i,t):
        '''Returns the total PAA for the right telecope'''
        calc_ang=LA.angle(self.v_r_func_tot(i,t),-self.u_r_func_tot(i,t))
        return calc_ang

    def calc_PAA_rin(self,i,t):
        '''Returns the inplane PAA for the right telecope'''
        calc_ang=LA.ang_in_out(self.v_r_func_tot(i,t),-self.u_r_func_tot(i,t),self.n_func(i,t),self.r_func(i,t),give='in')
        return calc_ang

    def calc_PAA_rout(self,i,t):
        '''Returns the outplane PAA for the right telecope'''
        calc_ang=LA.ang_in_out(self.v_r_func_tot(i,t),-self.u_r_func_tot(i,t),self.n_func(i,t),self.r_func(i,t),give='out')

        return calc_ang

    # Velocity
    def velocity_abs_calc(self,i_select,t,hstep):
        '''Returns the velocity vector of a spacecraft'''
        v = (np.array(self.putp(i_select,np.float64(t+hstep)))-np.array(self.putp(i_select,t)))/hstep
        return v

    def velocity_abs(self):
        '''Returns the velocity vector in a function'''
        hstep = np.float64(self.stat.hstep)
        #print(self.stat.hstep)
        #hstep = np.float128(hstep)
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

#    def high_precision(self,p):
#        '''Returns a high precision fit by the use of fourier components'''
#        Y = p
#        for i in range(0,len(p)):
#            y_p=[]
#            for j in range(0,len(p[i][0])):
#                y = p[i,:,j]
#                y_inv = scipy.fftpack.ifft(y)
#                y_new = scipy.fftpack.fft(y_inv)
#                Y[i,:,j] = np.real(y_new)
#        return Y

    



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
        print('')
        print('Importing Orbit')
        tic=time.clock()
        Orbit=orbit.ORBIT(input_param=self.stat.__dict__)
        try:
            print(str(Orbit.linecount)+' datapoints')
        except AttributeError:
            print('Obtained orbital function')
        self.orbit = Orbit
        self.LISA = self.orbit.LISA
        self.putp = self.orbit.putp
        print('Done in '+str(time.clock()-tic))
        self.SC = range(1,4) 
        #self.COM_func = const.COM_func(self)

        # Calculations
        v_l_func_tot=[]
        v_r_func_tot=[]
        u_l_func_tot=[]
        u_r_func_tot=[]
        #v_l0test_func_tot=[]
        #v_r0test_func_tot=[]
        #u_l0test_func_tot=[]
        #u_r0test_func_tot=[]
        L_sl_func_tot=[]
        L_sr_func_tot=[]
        L_rl_func_tot=[]
        L_rr_func_tot=[]
        v_l_stat_func_tot=[]
        v_r_stat_func_tot=[]
        
        #--- Obtaining Velocity
        self.velocity_abs()
        self.velocity_func()

        for i in range(1,4):
            [[v_l_func,v_r_func,u_l_func,u_r_func],[L_sl_func,L_sr_func,L_rl_func,L_rr_func],[v_l0_func,v_r0_func,u_l0_func,u_r0_func]] = self.send_func(i)

            v_l_func_tot.append(v_l_func)
            v_r_func_tot.append(v_r_func)
            u_l_func_tot.append(u_l_func)
            u_r_func_tot.append(u_r_func)
            #v_l0test_func_tot.append(v_l0_func)
            #v_r0test_func_tot.append(v_r0_func)
            #u_l0test_func_tot.append(u_l0_func)
            #u_r0test_func_tot.append(u_r0_func)
            
            L_sl_func_tot.append(L_sl_func)
            L_sr_func_tot.append(L_sr_func)
            L_rl_func_tot.append(L_rl_func)
            L_rr_func_tot.append(L_rr_func)
            
            #v_l_stat_func_tot.append(const.get_armvec_func(self,i_self,'l'))
            #v_r_stat_func_tot.append(const.get_armvec_func(self,i_self,'r'))
        
        self.v_l_func_tot = self.func_over_sc(v_l_func_tot)
        self.v_r_func_tot = self.func_over_sc(v_r_func_tot)
        self.u_l_func_tot = self.func_over_sc(u_l_func_tot)
        self.u_r_func_tot = self.func_over_sc(u_r_func_tot)
        #self.v_l0test_func_tot = const.func_over_sc(v_l0test_func_tot)
        #self.v_r0test_func_tot = const.func_over_sc(v_r0test_func_tot)
        #self.u_l0test_func_tot = const.func_over_sc(u_l0test_func_tot)
        #self.u_r0test_func_tot = const.func_over_sc(u_r0test_func_tot)

        self.L_sl_func_tot = self.func_over_sc(L_sl_func_tot)
        self.L_sr_func_tot = self.func_over_sc(L_sr_func_tot)
        self.L_rl_func_tot = self.func_over_sc(L_rl_func_tot)
        self.L_rr_func_tot = self.func_over_sc(L_rr_func_tot)
        
        #self.v_l_stat_func_tot = const.func_over_sc(v_l_stat_func_tot)
        #self.v_r_stat_func_tot = const.func_over_sc(v_r_stat_func_tot)
        self.v_l_stat_func_tot = lambda i,t: self.get_armvec_func(i,t,'l')
        self.v_r_stat_func_tot = lambda i,t: self.get_armvec_func(i,t,'r')

        self.n_func = lambda i,t: LA.unit(np.cross(self.v_l_stat_func_tot(i,t),self.v_r_stat_func_tot(i,t)))
        self.r_func = lambda i,t: (self.v_l_stat_func_tot(i,t)+self.v_r_stat_func_tot(i,t))/2.0
        #self.r_func = lambda i,t: const.r_calc(self.v_l_stat_func_tot(i,t),self.v_r_stat_func_tot(i,t),i)
        #self.pos_func = const.func_over_sc(pos_func)

        self.v_l_in_func_tot = lambda i,t: LA.inplane(self.v_l_func_tot(i,t),self.n_func(i,t))
        self.v_r_in_func_tot = lambda i,t: LA.inplane(self.v_r_func_tot(i,t),self.n_func(i,t))
        self.u_l_in_func_tot = lambda i,t: LA.inplane(self.u_l_func_tot(i,t),self.n_func(i,t))
        self.u_r_in_func_tot = lambda i,t: LA.inplane(self.u_r_func_tot(i,t),self.n_func(i,t))
        self.v_l_out_func_tot = lambda i,t: LA.outplane(self.v_l_func_tot(i,t),self.n_func(i,t))
        self.v_r_out_func_tot = lambda i,t: LA.outplane(self.v_r_func_tot(i,t),self.n_func(i,t))
        self.u_l_out_func_tot = lambda i,t: LA.outplane(self.u_l_func_tot(i,t),self.n_func(i,t))
        self.u_r_out_func_tot = lambda i,t: LA.outplane(self.u_r_func_tot(i,t),self.n_func(i,t))

        #--- Obtaining PAA --- 
        selections=['l_in','l_out','r_in','r_out']
        PAA_func_val={}
        PAA_func_val[selections[0]] = lambda i,t: self.calc_PAA_lin(i,t)
        PAA_func_val[selections[1]] = lambda i,t: self.calc_PAA_lout(i,t)
        PAA_func_val[selections[2]] = lambda i,t: self.calc_PAA_rin(i,t)
        PAA_func_val[selections[3]] = lambda i,t: self.calc_PAA_rout(i,t)
        PAA_func_val['l_tot'] = lambda i,t: const.calc_PAA_ltot(i,t)
        PAA_func_val['r_tot'] = lambda i,t: const.calc_PAA_rtot(i,t)

        self.PAA_func = PAA_func_val 
       
        self.ang_breathing_din = lambda i, time: LA.angle(self.v_l_func_tot(i,time),self.v_r_func_tot(i,time))
        self.ang_breathing_in = lambda i, time: LA.angle(self.u_l_func_tot(i,time),self.u_r_func_tot(i,time))
        self.ang_breathing_stat = lambda i, time: LA.angle(self.v_l_stat_func_tot(i,time),self.v_r_stat_func_tot(i,time))
        
        self.ang_in_l = lambda i,t: LA.ang_in(self.v_l_func_tot(i,t),self.n_func(i,t),self.r_func(i,t))
        self.ang_in_r = lambda i,t: LA.ang_in(self.v_r_func_tot(i,t),self.n_func(i,t),self.r_func(i,t))
        self.ang_out_l = lambda i,t: LA.ang_out(self.v_l_func_tot(i,t),self.n_func(i,t))
        self.ang_out_r = lambda i,t: LA.ang_out(self.v_r_func_tot(i,t),self.n_func(i,t))
        
        try:
            self.t_all
        except AttributeError:
            self.t_all = self.orbit.t
        if 'int' in str(type(self.stat.length_calc)):
            if self.t_all[-1]>(self.stat.length_calc+1)*self.param.day2sec:
                loc = calc.get_nearest_smaller_value(self.orbit.t,self.stat.length_calc*self.param.day2sec)
                self.t_all = self.orbit.t[0:loc+1]


        return self

    def solve_L_PAA(self,t,pos_self,pos_left,pos_right,select='sl',calc_method='Waluschka',i=False):
        '''Calculate the photon traveling time along one of the six laserlinks'''
        t_guess = np.linalg.norm(np.array(self.putp(1,0)) - np.array(self.putp(2,0)))/c

        if select=='sl' or select=='rl':
            s1 = lambda x: pos_left(x)
        elif select=='sr' or select=='rr':
            s1 = lambda x: pos_right(x)

        s2 = lambda x: pos_self(x)
        x_0 = t
        if select=='sl' or select=='sr':
            if calc_method=='Abram':
                s3 = lambda dt: s1(x_0+dt) - s2(x_0)
            elif calc_method=='Waluschka':
                s3 = lambda dt: s1(x_0+dt) - s2(x_0+dt)
        elif select=='rl' or select=='rr':
            if calc_method=='Abram':
                s3 = lambda dt: -s1(x_0-dt) + s2(x_0)
            elif calc_method=='Waluschka':
                s3 = lambda dt: -s1(x_0-dt) + s2(x_0-dt)
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
        calc_method = self.stat.calc_method
        selections = ['sl','sr','rl','rr']

        L_sl_func =  lambda time: self.solve_L_PAA(time,pos_self,pos_left,pos_right,select=selections[0],calc_method=calc_method,i=i)
        L_sr_func =  lambda time: self.solve_L_PAA(time,pos_self,pos_left,pos_right,select=selections[1],calc_method=calc_method,i=i)
        L_rl_func =  lambda time: self.solve_L_PAA(time,pos_self,pos_left,pos_right,select=selections[2],calc_method=calc_method,i=i)
        L_rr_func =  lambda time: self.solve_L_PAA(time,pos_self,pos_left,pos_right,select=selections[3],calc_method=calc_method,i=i)

        return [L_sl_func,L_sr_func,L_rl_func,L_rr_func]

    def relativistic_aberrations(self,i,t,u):
        '''Adjust vecor u by adding the angle caused by aberration'''
        relativistic=self.stat.relativistic

        if self.stat.calc_method=='Abram':
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

        elif self.stat.calc_method=='Waluschka':

            return u



#    def r_calc(self,v_l,v_r,i,m=[2,2,2]):
#        '''Returns the vector r pointing from a spacecraft towards the COMof the constellation'''
#        [i_OBJ,i_left,i_right] = self.i_slr(i)
#        r =  (v_l*m[i_left-1]+v_r*m[i_right-1])/(m[i_left-1]+m[i_right-1])
#
#        return r

    def send_func(self,i):
        '''Uses previous defined functions to return the vecors L, u, v, r and n'''
        calc_method=self.stat.calc_method
        [i_self,i_left,i_right] = utils.const.i_slr(i)

        pos_left = lambda time: self.putp(i_left,time)
        pos_self = lambda time: self.putp(i_self,time)
        pos_right = lambda time: self.putp(i_right,time)

        [L_sl,L_sr,L_rl,L_rr] = self.L_PAA(pos_self,pos_left,pos_right,i=i_self)

        if calc_method=='Abram':
            #Abram2018
            v_send_l0 = lambda t: pos_left(t+L_sl(t)) - pos_self(t)
            v_send_r0 = lambda t: pos_right(t+L_sr(t)) - pos_self(t)
            v_rec_l0 = lambda t: pos_self(t) - pos_left(t - L_rl(t))
            v_rec_r0 = lambda t: pos_self(t) - pos_right(t - L_rr(t))

        elif calc_method=='Waluschka':
            #Waluschka2003
            v_send_l0 = lambda t: pos_left(t+L_sl(t)) - pos_self(t+L_sl(t))
            v_send_r0 = lambda t: pos_right(t+L_sr(t)) - pos_self(t+L_sr(t))
            v_rec_l0 = lambda t: pos_self(t-L_rl(t)) - pos_left(t - L_rl(t))
            v_rec_r0 = lambda t: pos_self(t-L_rr(t)) - pos_right(t - L_rr(t))

        if self.stat.aberration==False:
            v_send_l = v_send_l0
            v_send_r = v_send_r0
            v_rec_l = v_rec_l0
            v_rec_r = v_rec_r0
        elif self.stat.aberration==True:
            v_send_l = lambda t: self.relativistic_aberrations(i,t,v_send_l0(t))
            v_send_r = lambda t: self.relativistic_aberrations(i,t,v_send_r0(t))
            v_rec_l = lambda t: self.relativistic_aberrations(i,t,v_rec_l0(t))
            v_rec_r = lambda t: self.relativistic_aberrations(i,t,v_rec_r0(t))

        return [[v_send_l,v_send_r,v_rec_l,v_rec_r],[L_sl,L_sr,L_rl,L_rr],[v_send_l0,v_send_r0,v_rec_l0,v_rec_r0]]
