from pointLISA import *   
import LA
import pointLISA.methods as methods
from synthlisa import *
import numpy as np
import matplotlib.pyplot as plt
import os
from fractions import Fraction
import math
import datetime
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
import scipy.optimize
from sympy import *
from imports import *
from parameters import *

# This class contains different helper functions

class Object(object):
    '''Creates a (empty) class'''
    pass

def nominal_arm(OBJ,i,t):
    '''Returns a functions of the normalized OBJ.orbit.L'''
    L_vec=[]
    t_vec=OBJ.orbit.t
    for j in range(0,len(t_vec)):
        L_vec.append(np.linalg.norm(orbit.L[i-1][j]))

    f=interp1d(t_vec,L_vec,bounds_error)

    return f(t)

def LISA_obj(OBJ,type_select='cache'):
    '''Creates the attribute LISA for OBJ which is a synthLISA object of a pre-setted type'''
    if 'function' in str(type(type_select)):
        lisa = Object()
        lisa.putp = type_select
        OBJ.LISA = lisa

    else:
        func_nominal_arm = lambda i,time: nominal_arm(OBJ,i,time)
        lisa=OBJ.orbit.lisa_obj
        lisa_orb=PyLISA(lisa,func_nominal_arm)
        lisa_cache=CacheLISA(lisa_orb)

        OBJ.t_all = OBJ.orbit.t
        if type(type_select)==str:
            if type_select=='cache':
                OBJ.LISA = lisa_cache
            elif type_select=='Py':
                OBJ.LISA = lisa_orb
            elif type_select=='other':
                OBJ.LISA = lisa
        else:
            OBJ.LISA=type_select
    return 0

def i_slr(i,side='all'):
    '''Obtaining the correct spacecraft numbers'''

    i_OBJ = i
    i_left = (i+1)%3
    i_right = (i+2)%3

    i_ret = [i_OBJ,i_left,i_right]
    for j in range(0,len(i_ret)):
        if i_ret[j]==0:
            i_ret[j]=3

    if side=='all':
        return i_ret
    elif side=='l':
        return [i_ret[0],i_ret[1]]
    elif side=='r':
        return [i_ret[0],i_ret[2]]

#def get_link(i,side):
#    if side=='l':
#        link=(i+2)%3
#    if side=='r':
#        link=(i+1)%3
#
#    if link==0:
#        link=3
#
#    return link

def get_armvec_func(OBJ,i,side):
    '''Obtains the functions of the distance vectors between two spacecrafts'''
    [i_OBJ,i_next] = i_slr(i,side=side)
    arm_vec = lambda time: np.array(OBJ.putp(i_next,time)) - np.array(OBJ.putp(i_OBJ,time))

    return arm_vec

def func_pos(OBJ,i):
    '''Generate functions of the positions''' 
    if OBJ.test_COM_effect == False:
        L = lambda time: np.array(OBJ.putp(i,time))
    if OBJ.test_COM_effect == True:
        L = lambda time: np.array(OBJ.putp(i,time) - OBJ.COM_func(time))
    return L

def COM_func(OBJ):
    '''This function obtaines the coordinate function of the center of mass of the LISA constellation'''
    COM = lambda time: (OBJ.putp(1,time)+OBJ.putp(2,time)+OBJ.putp(3,time))/3.0 
    return COM

def solve_L_PAA(OBJ,t,pos_OBJ,pos_left,pos_right,select='sl',calc_method='Waluschka',i=False):
    '''Calculate the photon traveling time along one of the six laserlinks'''
    if OBJ.LISA==False:
        t_guess = np.linalg.norm(OBJ.orbit.p[0][0,:] - OBJ.orbit.p[1][0,:])/c
    else:
        t_guess = np.linalg.norm(np.array(OBJ.putp(1,0)) - np.array(OBJ.putp(2,0)))/c
    
    if OBJ.test_COM_effect==False:
        if select=='sl' or select=='rl':
            s1 = lambda x: pos_left(x)
        elif select=='sr' or select=='rr':
            s1 = lambda x: pos_right(x)

        s2 = lambda x: pos_OBJ(x)
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

        res = scipy.optimize.brentq(s5,0,t_guess*4)
    
    elif OBJ.test_COM_effect==True:
        if select=='sl' or select=='rl':
            s1 = lambda x: pos_left(x)
        elif select=='sr' or select=='rr':
            s1 = lambda x: pos_right(x)

        s2 = lambda x: pos_OBJ(x)
        x_0 = t
        if select=='sl' or select=='sr':
            if calc_method=='Abram':
                s3 = lambda dt: s1(x_0+dt) - s2(x_0)
                com = lambda dt: OBJ.COM_func(x_0+dt) - OBJ.COM_func(x_0)
            elif calc_method=='Waluschka':
                s3 = lambda dt: s1(x_0+dt) - s2(x_0+dt)
                com = lambda dt: OBJ.COM_func(x_0+dt) - OBJ.COM_func(x_0+dt)
        elif select=='rl' or select=='rr':
            if calc_method=='Abram':
                s3 = lambda dt: -s1(x_0-dt) + s2(x_0)
                com = lambda dt: -OBJ.COM_func(x_0-dt) + OBJ.COM_func(x_0)
            elif calc_method=='Waluschka':
                s3 = lambda dt: -s1(x_0-dt) + s2(x_0-dt)
                com = lambda dt: -OBJ.COM_func(x_0-dt) + OBJ.COM_func(x_0-dt)
        s4 = lambda dt: np.linalg.norm(s3(dt)-com(dt))
        s5 = lambda dt: s4(dt) - c*dt
        
        res = scipy.optimize.brentq(s5,0,t_guess*4)

    return res


def L_PAA(OBJ,pos_OBJ,pos_left,pos_right,calc_method='Walushka',i=False):
    '''Obtain time of flight of beam between spacecrafts'''

    selections = ['sl','sr','rl','rr']

    L_sl_func =  lambda time: solve_L_PAA(OBJ,time,pos_OBJ,pos_left,pos_right,select=selections[0],calc_method=calc_method,i=i)
    L_sr_func =  lambda time: solve_L_PAA(OBJ,time,pos_OBJ,pos_left,pos_right,select=selections[1],calc_method=calc_method,i=i)
    L_rl_func =  lambda time: solve_L_PAA(OBJ,time,pos_OBJ,pos_left,pos_right,select=selections[2],calc_method=calc_method,i=i)
    L_rr_func =  lambda time: solve_L_PAA(OBJ,time,pos_OBJ,pos_left,pos_right,select=selections[3],calc_method=calc_method,i=i)

    return [L_sl_func,L_sr_func,L_rl_func,L_rr_func]

def r_calc(v_l,v_r,i,m=[2,2,2]):
    '''Returns the vector r pointing from a spacecraft towards the COMof the constellation'''
    [i_OBJ,i_left,i_right] = i_slr(i)
    r =  (v_l*m[i_left-1]+v_r*m[i_right-1])/(m[i_left-1]+m[i_right-1])

    return r

def func_over_sc(func_tot):
    '''Makes from a list of funcions a function (wih two variables'''
    f = lambda i,t: func_tot[i-1](t)

    return f

def send_func(OBJ,i,calc_method='Waluschka'):
    '''Uses previous defined functions to return the vecors L, u, v, r and n'''
    [i_OBJ,i_left,i_right] = i_slr(i)

    pos_left = func_pos(OBJ,i_left)
    pos_OBJ = func_pos(OBJ,i_OBJ)
    pos_right = func_pos(OBJ,i_right)

    if OBJ.delay==True:
        [L_sl,L_sr,L_rl,L_rr] = L_PAA(OBJ,pos_OBJ,pos_left,pos_right,calc_method=calc_method,i=i_OBJ)
    elif OBJ.delay=='not ahead':
        L_sl = lambda t: np.linalg.norm(pos_left(t) - pos_OBJ(t))/c
        L_sr = lambda t: np.linalg.norm(pos_right(t) - pos_OBJ(t))/c
        L_rl=L_sl
        L_rr=L_sr

    elif OBJ.delay=='constant':
        L_sl = lambda t: OBJ.armlength/c #...adjust
        L_sr = lambda t: OBJ.armlength/c
        L_rl=L_sl
        L_rr=L_sr


    elif OBJ.delay==False:
        L_sl = lambda t: 0
        L_sr = lambda t: 0
        L_rl=L_sl
        L_rr=L_sr
    
    if OBJ.test_COM_effect==False:
        if calc_method=='Abram':
            #Abram2018
            v_send_l0 = lambda t: pos_left(t+L_sl(t)) - pos_OBJ(t)
            v_send_r0 = lambda t: pos_right(t+L_sr(t)) - pos_OBJ(t)
            v_rec_l0 = lambda t: pos_OBJ(t) - pos_left(t - L_rl(t))
            v_rec_r0 = lambda t: pos_OBJ(t) - pos_right(t - L_rr(t))

        elif calc_method=='Waluschka':
            #Waluschka2003
            v_send_l0 = lambda t: pos_left(t+L_sl(t)) - pos_OBJ(t+L_sl(t))
            v_send_r0 = lambda t: pos_right(t+L_sr(t)) - pos_OBJ(t+L_sr(t))
            v_rec_l0 = lambda t: pos_OBJ(t-L_rl(t)) - pos_left(t - L_rl(t))
            v_rec_r0 = lambda t: pos_OBJ(t-L_rr(t)) - pos_right(t - L_rr(t))

    elif OBJ.test_COM_effect==True:
        if calc_method=='Abram':
            #Abram2018
            v_send_l0 = lambda t: pos_left(t+L_sl(t)) - pos_OBJ(t) - (OBJ.COM_func(t+L_sl(t)) - OBJ.COM_func(t))
            v_send_r0 = lambda t: pos_right(t+L_sr(t)) - pos_OBJ(t) - (OBJ.COM_func(t+L_sr(t)) - OBJ.COM_func(t))
            v_rec_l0 = lambda t: pos_OBJ(t) - pos_left(t - L_rl(t)) - (OBJ.COM_func(t) - OBJ.COM_func(t-L_rl(t)))
            v_rec_r0 = lambda t: pos_OBJ(t) - pos_right(t - L_rr(t)) - (OBJ.COM_func(t) - OBJ.COM_func(t-L_rr(t)))

        elif calc_method=='Waluschka':
            #Waluschka2003
            v_send_l0 = lambda t: pos_left(t+L_sl(t)) - pos_OBJ(t+L_sl(t))
            v_send_r0 = lambda t: pos_right(t+L_sr(t)) - pos_OBJ(t+L_sr(t))
            v_rec_l0 = lambda t: pos_OBJ(t-L_rl(t)) - pos_left(t - L_rl(t))
            v_rec_r0 = lambda t: pos_OBJ(t-L_rr(t)) - pos_right(t - L_rr(t))

    if OBJ.aberration==False:
        v_send_l = v_send_l0
        v_send_r = v_send_r0
        v_rec_l = v_rec_l0
        v_rec_r = v_rec_r0
    elif OBJ.aberration==True:
        v_send_l = lambda t: relativistic_aberrations(OBJ,i,t,v_send_l0(t),relativistic=OBJ.relativistic)
        v_send_r = lambda t: relativistic_aberrations(OBJ,i,t,v_send_r0(t),relativistic=OBJ.relativistic)
        v_rec_l = lambda t: relativistic_aberrations(OBJ,i,t,v_rec_l0(t),relativistic=OBJ.relativistic)
        v_rec_r = lambda t: relativistic_aberrations(OBJ,i,t,v_rec_r0(t),relativistic=OBJ.relativistic)

    return [[v_send_l,v_send_r,v_rec_l,v_rec_r],[L_sl,L_sr,L_rl,L_rr],[v_send_l0,v_send_r0,v_rec_l0,v_rec_r0]]

def relativistic_aberrations(OBJ,i,t,u,relativistic=True):
    '''Adjust vecor u by adding the angle caused by aberration'''
    if OBJ.calc_method=='Abram':
        if relativistic==True:
            V = -OBJ.vel.abs(i,t)
            V_mag = np.linalg.norm(V)
            u_mag = np.linalg.norm(u)
            c_vec = LA.unit(u)*c

            velo = V
            coor = methods.coor_SC(OBJ,i,t)
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

    #        ux = (np.dot(c_vec,V)/(V_mag))*LA.unit(V)
    #        x = LA.unit(ux)
    #        uy = c_vec-ux
    #        y = LA.unit(uy)
    #        
    #        ux_mag = np.linalg.norm(ux)
    #        uy_mag = np.linalg.norm(uy)
    #
    #        den = (1+((ux_mag*V_mag)/(c**2)))
    #        ux_ac = (ux_mag+V_mag)/den
    #        gamma = 1.0/((1-((V_mag/c)**2))**0.5)
    #        uy_ac = uy_mag/(gamma*den)
    #        
    #        u_new = LA.unit(ux_ac*x+uy_ac*y)*u_mag
    #        #u_new = (ux_ac*x+uy_ac*y)*(u_mag/c)

        elif relativistic==False:
            V = -OBJ.vel.abs(i,t)
            u_mag = np.linalg.norm(u)
            c_vec = LA.unit(u)*c

            u_new = LA.unit(c_vec+V)*u_mag
        else:
            print('Error')

        return u_new

    elif OBJ.calc_method=='Waluschka':

        return u



#def relativistic_aberrations2(OBJ,i,t,tdel,side,relativistic=True): #not used
#    '''Adjust vecor u by adding the angle caused by aberration'''
#    [i_self,i_left,i_right] = i_slr(i)
#    if OBJ.calc_method=='Abram':
#        tdel0=0
#    elif OBJ.calc_method=='Waluschka':
#        tdel0 = tdel
#
#    if side=='l':
#        u_not_ab = np.array(OBJ.putp(i_self,t-tdel0)) - np.array(OBJ.putp(i_left,t-tdel))
#        u_ab = np.linalg.norm(u_not_ab)*(LA.unit(LA.unit(u_not_ab)*OBJ.c+(OBJ.vel.abs(i_self,t-tdel0) - OBJ.vel.abs(i_left,t-tdel))))
#
#    elif side=='r':
#        u_not_ab = np.array(OBJ.putp(i_self,t-tdel0)) - np.array(OBJ.putp(i_right,t-tdel))
#        u_ab = np.linalg.norm(u_not_ab)*(LA.unit(LA.unit(u_not_ab)*OBJ.c+(OBJ.vel.abs(i_self,t-tdel0) - OBJ.vel.abs(i_right,t-tdel))))
# 
#    if relativistic==False:
#        return u_ab
#    
#    elif relativistic==True:
#        coor = methods.coor_SC(OBJ,i_self,t-tdel0)
#        if side=='l':
#            velo = (OBJ.vel.abs(i_self,t-tdel0) - OBJ.vel.abs(i_left,t-tdel))
#        elif side=='r':
#            velo = (OBJ.vel.abs(i_self,t-tdel0) - OBJ.vel.abs(i_right,t-tdel))
#
#        c_vec = LA.unit(u_not_ab)*c
#
#        r = coor[0]
#        x_prime = LA.unit(velo)
#        n_prime = LA.unit(np.cross(velo,r))
#        r_prime = LA.unit(np.cross(n_prime,x_prime))
#
#        coor_velo = np.array([r_prime,n_prime,x_prime])
#        c_velo = LA.matmul(coor_velo,c_vec)
#        v = np.linalg.norm(velo)
#        den = 1.0 - ((v/(c**2))*coor_velo[2])
#        num = ((1.0-((v**2)/(c**2)))**0.5)
#
#        ux_prime = (c_velo[2] - v)/den
#        ur_prime = (num*c_velo[0])/den
#        un_prime = (num*c_velo[1])/den
#        c_prime = ux_prime*x_prime + un_prime*n_prime +ur_prime*r_prime
#        u_new=LA.unit(c_prime)*np.linalg.norm(u_not_ab)
#   
#    return u_new

#PAA angles
def calc_PAA_ltot(OBJ,i,t):
    '''Returns the total PAA for the left telecope'''
    calc_ang=LA.angle(OBJ.v_l_func_tot(i,t),-OBJ.u_l_func_tot(i,t))
    return calc_ang

def calc_PAA_lin(OBJ,i,t):
    '''Returns the inplane PAA for the left telecope'''
    calc_ang=LA.ang_in_out(OBJ.v_l_func_tot(i,t),-OBJ.u_l_func_tot(i,t),OBJ.n_func(i,t),OBJ.r_func(i,t),give='in')
    return calc_ang

def calc_PAA_lout(OBJ,i,t):
    '''Returns the outplanr PAA for the left telecope'''
    calc_ang=LA.ang_in_out(OBJ.v_l_func_tot(i,t),-OBJ.u_l_func_tot(i,t),OBJ.n_func(i,t),OBJ.r_func(i,t),give='out')
    return calc_ang

def calc_PAA_rtot(OBJ,i,t):
    '''Returns the total PAA for the right telecope'''
    calc_ang=LA.angle(OBJ.v_r_func_tot(i,t),-OBJ.u_r_func_tot(i,t))
    return calc_ang

def calc_PAA_rin(OBJ,i,t):
    '''Returns the inplane PAA for the right telecope'''
    calc_ang=LA.ang_in_out(OBJ.v_r_func_tot(i,t),-OBJ.u_r_func_tot(i,t),OBJ.n_func(i,t),OBJ.r_func(i,t),give='in')
    return calc_ang

def calc_PAA_rout(OBJ,i,t):
    '''Returns the outplane PAA for the right telecope'''
    calc_ang=LA.ang_in_out(OBJ.v_r_func_tot(i,t),-OBJ.u_r_func_tot(i,t),OBJ.n_func(i,t),OBJ.r_func(i,t),give='out')
    
    return calc_ang

# Velocity
def velocity_abs_calc(OBJ,i_select,t,hstep):
    '''Returns the velocity vector of a spacecraft'''
    v = (np.array(OBJ.putp(i_select,np.float64(t+hstep)))-np.array(OBJ.putp(i_select,t)))/hstep
    return v


def velocity_abs(OBJ,hstep=1.0):
    '''Returns the velocity vector in a function'''
    hstep = np.float128(hstep)
    v_ret = lambda i,time: velocity_abs_calc(OBJ,i,time,hstep)
    try:
        OBJ.vel
    except AttributeError:
        OBJ.vel = Object()
    OBJ.vel.abs = v_ret
    
    return OBJ.vel.abs
    
def velocity_calc(OBJ,i,time,hstep,side,rs):
    '''Calculates the velocity components'''
    [i_OBJ,i_next] = i_slr(i,side=side)
    v_pos_l = OBJ.v_l_stat_func_tot(i_OBJ,time)
    v_pos_r = OBJ.v_r_stat_func_tot(i_OBJ,time)
    n = np.cross(v_pos_r,v_pos_l)
    n = LA.unit(n)
    if side=='l':
        v_pos = v_pos_l
    elif side=='r':
        v_pos = v_pos_r 

    pos_OBJ = np.array(OBJ.putp(i_OBJ,time))
    pos_next = np.array(OBJ.putp(i_next,time))
    pos_OBJ_h = np.array(OBJ.putp(i_OBJ,time+np.float64(hstep)))
    pos_next_h = np.array(OBJ.putp(i_next,time+np.float64(hstep)))
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

def velocity_func(OBJ,hstep=1.0):
    '''Returns functions of al velocity components'''
    hstep = np.float64(hstep)
    
    try:
        OBJ.vel
    except AttributeError:
        OBJ.vel = Object()
    OBJ.vel.l= lambda i,time: velocity_calc(OBJ,i,time,hstep,'l',0)
    OBJ.vel.in_l = lambda i,time: velocity_calc(OBJ,i,time,hstep,'l',1)
    OBJ.vel.out_l = lambda i,time: velocity_calc(OBJ,i,time,hstep,'l',2)
    OBJ.vel.arm_l = lambda i,time: velocity_calc(OBJ,i,time,hstep,'l',3)
    OBJ.vel.in_mag_l = lambda i,time: velocity_calc(OBJ,i,time,hstep,'l',4)
    OBJ.vel.out_mag_l = lambda i,time: velocity_calc(OBJ,i,time,hstep,'l',5)
    OBJ.vel.arm_mag_l = lambda i,time: velocity_calc(OBJ,i,time,hstep,'l',6)
    OBJ.vel.r= lambda i,time: velocity_calc(OBJ,i,time,hstep,'r',0)
    OBJ.vel.in_r = lambda i,time: velocity_calc(OBJ,i,time,hstep,'r',1)
    OBJ.vel.out_r = lambda i,time: velocity_calc(OBJ,i,time,hstep,'r',2)
    OBJ.vel.arm_r = lambda i,time: velocity_calc(OBJ,i,time,hstep,'r',3)
    OBJ.vel.in_mag_r = lambda i,time: velocity_calc(OBJ,i,time,hstep,'r',4)
    OBJ.vel.out_mag_r = lambda i,time: velocity_calc(OBJ,i,time,hstep,'r',5)
    OBJ.vel.arm_mag_r = lambda i,time: velocity_calc(OBJ,i,time,hstep,'r',6)
    
    return 0


def high_precision(p):
    '''Returns a high precision fit by the use of fourier components'''
    Y = p
    for i in range(0,len(p)):
        y_p=[]
        for j in range(0,len(p[i][0])):
            y = p[i,:,j]
            y_inv = scipy.fftpack.ifft(y)
            y_new = scipy.fftpack.fft(y_inv)
            Y[i,:,j] = np.real(y_new)
    return Y
