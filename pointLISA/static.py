##!/usr/bin/env python 
from pointLISA import *

# in the STAT class an ORBIT object is used which containes the imported or constructed orbital coordinates of the spacecrafts. In STAT different vectors are being constructed and calculaded (r, n, v, u, L) as well as angles (PAA, breathing)

class STAT():
    def __init__(self,input_param,**kwargs):
        self.stat = self.get_settings(input_param,kwargs)
        if 'function' not in str(type(self.stat.LISA_opt)):
            self.get_scale_and_time_units()

    def get_settings(self,input_param,kwargs):
        stat = utils.Object()

        for k in input_param.keys():
            setattr(stat,k,input_param[k])

        for key,value in kwargs.items():
            input_param[key] = value
            setattr(stat,key,value)

        return stat

    def get_scale_and_time_units(self):
        if self.stat.scale=='Default':
            print('Getting scale by filename:')
            a = self.stat.filename
            a1 = a.split('.')[0]
            a1 = a1.split('_')
            for k in range(0,len(a1)):
                if 'scale' == a1[k]:
                    self.stat.scale = float(a1[k+1]) 
            print(self.stat.scale)
        print('')
        
        if self.stat.timeunit=='Default':
            print('Getting timestep by filename:')
            a = self.stat.filename
            a1 = a.split('.')[0]
            a1 = a1.split('_')
            for k in range(0,len(a1)):
                if 'timestep' == a1[k]:
                    self.stat.timeunit = a1[k+1]
                    print(self.stat.timeunit)
            if self.stat.timeunit!='days' and self.stat.timeunit!='seconds':
                print('Could not obtain proper timestep')
        print('')
 

    def putp(self,i,t,mode='Default'):
        '''Returns the coordinates of spacecraft i from a linear interpolation or by the putp function in synthLISA'''
        if mode=='Default':
            mode=self.stat.putp_mode
        if 'function' in str(type(self.stat.LISA_opt)):
            return self.LISA.putp(i,t)
        else:
            return self.putp_fitted(i,t)


    def PAA_func(self):
        '''Obtains functions of vectors and angles (PAA, brething, n, r, u, v, L'''
        print('')
        print('Importing Orbit')
        tic=time.clock()
        Orbit=orbit.ORBIT(input_param=self.stat.__dict__)
        print(str(Orbit.linecount)+' datapoints')
        self.orbit = Orbit
        utils.LISA_obj(self)
        print('Done in '+str(time.clock()-tic))
        self.SC = range(1,4) 
        self.putp_fitted = calc.get_putp_fitted(self)
        self.COM_func = utils.COM_func(self)

        # Calculations
        v_l_func_tot=[]
        v_r_func_tot=[]
        u_l_func_tot=[]
        u_r_func_tot=[]
        v_l0test_func_tot=[]
        v_r0test_func_tot=[]
        u_l0test_func_tot=[]
        u_r0test_func_tot=[]
        L_sl_func_tot=[]
        L_sr_func_tot=[]
        L_rl_func_tot=[]
        L_rr_func_tot=[]
        v_l_stat_func_tot=[]
        v_r_stat_func_tot=[]
        pos_func=[]
        
        #--- Obtaining Velocity
        utils.velocity_abs(self)
        utils.velocity_func(self)

        for i in range(1,4):
            [[v_l_func,v_r_func,u_l_func,u_r_func],[L_sl_func,L_sr_func,L_rl_func,L_rr_func],[v_l0_func,v_r0_func,u_l0_func,u_r0_func]] = utils.send_func(self,i)

            v_l_func_tot.append(v_l_func)
            v_r_func_tot.append(v_r_func)
            u_l_func_tot.append(u_l_func)
            u_r_func_tot.append(u_r_func)
            v_l0test_func_tot.append(v_l0_func)
            v_r0test_func_tot.append(v_r0_func)
            u_l0test_func_tot.append(u_l0_func)
            u_r0test_func_tot.append(u_r0_func)
            
            L_sl_func_tot.append(L_sl_func)
            L_sr_func_tot.append(L_sr_func)
            L_rl_func_tot.append(L_rl_func)
            L_rr_func_tot.append(L_rr_func)
            
            [i_self,i_left,i_right] = utils.i_slr(i)
            v_l_stat_func_tot.append(utils.get_armvec_func(self,i_self,'l'))
            v_r_stat_func_tot.append(utils.get_armvec_func(self,i_self,'r'))
            pos_func.append(utils.func_pos(self,i))

        self.v_l_func_tot = utils.func_over_sc(v_l_func_tot)
        self.v_r_func_tot = utils.func_over_sc(v_r_func_tot)
        self.u_l_func_tot = utils.func_over_sc(u_l_func_tot)
        self.u_r_func_tot = utils.func_over_sc(u_r_func_tot)
        self.v_l0test_func_tot = utils.func_over_sc(v_l0test_func_tot)
        self.v_r0test_func_tot = utils.func_over_sc(v_r0test_func_tot)
        self.u_l0test_func_tot = utils.func_over_sc(u_l0test_func_tot)
        self.u_r0test_func_tot = utils.func_over_sc(u_r0test_func_tot)

        self.L_sl_func_tot = utils.func_over_sc(L_sl_func_tot)
        self.L_sr_func_tot = utils.func_over_sc(L_sr_func_tot)
        self.L_rl_func_tot = utils.func_over_sc(L_rl_func_tot)
        self.L_rr_func_tot = utils.func_over_sc(L_rr_func_tot)
        
        self.v_l_stat_func_tot = utils.func_over_sc(v_l_stat_func_tot)
        self.v_r_stat_func_tot = utils.func_over_sc(v_r_stat_func_tot)
        self.n_func = lambda i,t: LA.unit(np.cross(self.v_l_stat_func_tot(i,t),self.v_r_stat_func_tot(i,t)))
        self.r_func = lambda i,t: utils.r_calc(self.v_l_stat_func_tot(i,t),self.v_r_stat_func_tot(i,t),i)
        self.pos_func = utils.func_over_sc(pos_func)

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
        PAA_func_val[selections[0]] = lambda i,t: utils.calc_PAA_lin(self,i,t)
        PAA_func_val[selections[1]] = lambda i,t: utils.calc_PAA_lout(self,i,t)
        PAA_func_val[selections[2]] = lambda i,t: utils.calc_PAA_rin(self,i,t)
        PAA_func_val[selections[3]] = lambda i,t: utils.calc_PAA_rout(self,i,t)
        PAA_func_val['l_tot'] = lambda i,t: utils.calc_PAA_ltot(self,i,t)
        PAA_func_val['r_tot'] = lambda i,t: utils.calc_PAA_rtot(self,i,t)

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
            if self.t_all[-1]>(self.stat.length_calc+1)*day2sec:
                loc = calc.get_nearest_smaller_value(self.orbit.t,self.stat.length_calc*day2sec)
                self.t_all = self.orbit.t[0:loc+1]


        return self
