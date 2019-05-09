from imports import *
from pointLISA import *
import pointLISA
import LA
import numpy as np
import time
from orbit import ORBIT 
import utils
#year2sec=32536000
#day2sec=year2sec/365.25
#c=300000000

class STAT():
    def __init__(self,input_param,para,**kwargs):
        from imports import *
        for k in para:
            globals()[k] = para[k]
            setattr(self,k,para[k])
        for k in input_param:
            setattr(self,k,input_param[k])

        for key,value in kwargs.items():
            input_param[key] = value
            setattr(self,key,value)

        if self.scale=='Default':
            print('Getting scale by filename:')
            a = self.filename
            a1 = a.split('.')[0]
            a1 = a1.split('_')
            for k in range(0,len(a1)):
                if 'scale' == a1[k]:
                    self.scale = float(a1[k+1]) 
            print(self.scale)
        else:
            print(self.scale)
        input_param['scale']=self.scale
        print('')
        
        if self.timeunit=='Default':
            print('Getting timestep by filename:')
            a = self.filename
            a1 = a.split('.')[0]
            a1 = a1.split('_')
            for k in range(0,len(a1)):
                if 'timestep' == a1[k]:
                    self.timeunit = a1[k+1]
                    print(self.timeunit)
            if self.timeunit!='days' and self.timeunit!='seconds':
                print('Could not obtain proper timestep')
        else:
            print(self.timeunit)
        print('')
        input_param['timeunit'] = self.timeunit

        self.input_param = input_param

    def putp(self,i,t,mode='sampled'):
        if mode=='sampled':
            return self.putp_sampled(i,t)
        elif mode=='LISA':
            return self.LISA.putp(i,t)




    def PAA_func(self):
        print('')
        print('Importing Orbit')
        tic=time.clock()
        Orbit=ORBIT(input_param=self.input_param)
        print(str(Orbit.linecount)+' datapoints')
        self.orbit = Orbit
        utils.LISA_obj(self,type_select=self.LISA_opt)
        print('Done in '+str(time.clock()-tic))
        self.SC = range(1,4)
        
        self.putp_sampled = pointLISA.methods.get_putp_sampled(self)

        # Calculations
        #LA=utils.la()
        v_l_func_tot=[]
        v_r_func_tot=[]
        u_l_func_tot=[]
        u_r_func_tot=[]
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
        utils.velocity_func(self,hstep=self.hstep)
        utils.velocity_abs(self,hstep=self.hstep)

        for i in range(1,4):
            [[v_l_func,v_r_func,u_l_func,u_r_func],[L_sl_func,L_sr_func,L_rl_func,L_rr_func],[u_l0_func,u_r0_func]] = utils.send_func(self,i,calc_method = self.calc_method)

            v_l_func_tot.append(v_l_func)
            v_r_func_tot.append(v_r_func)
            u_l_func_tot.append(u_l_func)
            u_r_func_tot.append(u_r_func)
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
        print('Aberration: '+str(self.aberration))
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
        self.ang_breathing_stat = lambda i, time: LA.angle(self.v_l_stat_func_tot(i,time),self.v_r_stat_func_tot(i,time))
        
        self.ang_in_l = lambda i,t: LA.ang_in(self.v_l_func_tot(i,t),self.n_func(i,t),self.r_func(i,t))
        self.ang_in_r = lambda i,t: LA.ang_in(self.v_r_func_tot(i,t),self.n_func(i,t),self.r_func(i,t))
        self.ang_out_l = lambda i,t: LA.ang_out(self.v_l_func_tot(i,t),self.n_func(i,t))
        self.ang_out_r = lambda i,t: LA.ang_out(self.v_r_func_tot(i,t),self.n_func(i,t))

        return self


