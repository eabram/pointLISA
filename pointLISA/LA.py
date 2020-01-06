#from synthlisa import *
#import numpy as np
#import matplotlib.pyplot as plt
#import os
#from fractions import Fraction
#import math
#import datetime
#from scipy.interpolate import interp1d
#from scipy.interpolate import RegularGridInterpolator
#import scipy.optimize
#from sympy import *
from imports import *
import numpy as np
 
#year2sec=32536000
#day2sec=year2sec/365.25
#c=300000000

def norm(v):
    return np.linalg.norm(v)

def unit(v):
    if norm(v)==0:
        return v #...adjust
        raise ValueError
    else:
        return v/norm(v)

def angle(v1,v2,dot=False):
    norm_v1 = norm(v1)
    norm_v2 = norm(v2)
    if norm_v1!=0 and norm_v2!=0:
        if dot==False:
            sin = norm(np.cross(v1,v2)/(norm_v1*norm_v2))
            return np.arcsin(sin)
        elif dot == True:
            cos = np.dot(v1,v2)/(norm_v1*norm_v2)
            return np.sign(np.dot(v1,v2))*np.arccos(cos)
    else:
        #print('norm v1: '+str(norm_v1))
        #print('norm v2: '+str(norm_v2))

        return np.nan

#def aberration(thetas,v):
#    return (np.cos(thetas) - v/c)/(1-(v/c)*np.cos(thetas))

def inplane(v,n):
    inplane_calc = v - (np.dot(v,n)/(np.linalg.norm(n)**2))*n
    return inplane_calc

def outplane(v,n):
    outplane_calc = (np.dot(v,n)/(np.linalg.norm(n)**2))*n
    return outplane_calc

def ang_out(v,n):
    sign = np.sign(np.dot(outplane(v,n),n))
    return sign*angle(inplane(v,n),v)

def ang_in(v,n,r):
    #ang_out_calc = OBJ.ang_out(v,n)
    inplane_calc = inplane(v,n)
    ang_in_calc = angle(inplane_calc,r)

    return ang_in_calc

def ang_in_dot(v,v_stat,n,r):
    inplane_calc = inplane(v,n)
    costheta_beam = np.dot(inplane_calc,r)/(norm(inplane_calc)*norm(r))
    inplane_calc = inplane(v_stat,n)
    costheta_sc = np.dot(inplane_calc,r)/(norm(inplane_calc)*norm(r))

    return np.abs(np.arccos(costheta_beam)) - np.abs(np.arccos(costheta_sc))


def ang_in_direct(v,v_stat,n,r):
    inplane_calc = inplane(v,n)
    inplane_stat = inplane(v_stat,n)
    ang_out_calc = angle(inplane_calc,inplane_stat)
    #ang1 = OBJ.angle(inplane_calc,r)
    #ang2 = OBJ.angle(inplane_stat,r)
    #sign = 1#np.sign(ang1 - ang2)

    return ang_out_calc#ang1-ang2

def print_component(v,v_in,v_out,v_arm):
    n = norm(v)
    n_in = norm(v_in)
    n_out = norm(v_out)
    n_arm = norm(v_arm)

    print(n_in/n)
    print((n_out**2+n_in**2+n_arm**2)/n**2)
    print('')

    return 0

def ang_in_out(v1,v2,n,r,give='all'):
    n = unit(n)
    v1_out = (np.dot(v1,n)*n)/(norm(n)**2)
    #v1_arm = (np.dot(v1,v_stat)*v_stat)/(OBJ.norm(v_stat)**2)
    #v1_in = v1 - v1_out - v1_arm

    v2_out = (np.dot(v2,n)*n)/(norm(n)**2)
    #v2_arm = (np.dot(v2,v_stat)*v_stat)/(OBJ.norm(v_stat)**2)
    #v2_in = v2 - v2_out - v2_arm

    ang_out_1 = np.arcsin(norm(v1_out)/norm(v1))
    ang_out_1 = ang_out_1 * np.sign(np.dot(v1_out,n))
    ang_out_2 = np.arcsin(norm(v2_out)/norm(v2))
    ang_out_2 = ang_out_2 * np.sign(np.dot(v2_out,n))
    #ang_out = ang_out_1 - ang_out_2

    #q = np.cross(n,v_stat)
    #q = OBJ.unit(q)

    v1_in = v1 - v1_out
    v2_in = v2 - v2_out

    ang_in_1 = angle(v1_in,r)
    ang_in_2 = angle(v2_in,r)
    ang_in = ang_in_1 - ang_in_2
    ang_out = ang_out_1 - ang_out_2
    #ang_in = np.arcsin(np.linalg.norm(np.cross(v1_in,v2_in))/(np.linalg.norm(v1_in)*np.linalg.norm(v2_in)))


    #v1_in_calc = v1 - v1_out
    #v2_in_calc = v2 - v2_out
    #ang_in_1 = np.arcsin(OBJ.norm(np.cross(v1_in_calc,r))/(OBJ.norm(v1_in_calc)*OBJ.norm(r)))
    #ang_in_2 = np.arcsin(OBJ.norm(np.cross(v2_in_calc,r))/(OBJ.norm(v2_in_calc)*OBJ.norm(r)))


    #ang_in = OBJ.norm(np.cross(v1_in_calc,v2_in_calc))/(OBJ.norm(v1_in_calc)*OBJ.norm(v2_in_calc))
    #ang_in_1 = OBJ.norm(np.cross(v1_in_calc,v_stat))/(OBJ.norm(v1_in_calc)*OBJ.norm(v_stat))
    #ang_in_2 = OBJ.norm(np.cross(v2_in_calc,v_stat))/(OBJ.norm(v2_in_calc)*OBJ.norm(v_stat))

    #ang_in = abs(ang_in_1)+abs(ang_in_2)


    if give=='all':
        return [ang_in,ang_out]
    elif give=='in':
        return ang_in
    elif give=='out':
        return ang_out


def rotate(v,n,ang,mag=False):
    R = np.empty((3,3))
    c=np.cos(ang)
    s = np.sin(ang)
    [x,y,z]=n
    R[0,0] = c+(x**2)*(1-c)
    R[0,1] = x*y*(1-c) - z*s
    R[0,2] = x*z*(1-c)+y*s
    R[1,0] = y*x*(1-c) + z*s
    R[1,1] = c+(y**2)*(1-c)
    R[1,2] = y*z*(1-c)-x*s
    R[2,0] = z*x*(1-c) - y*s
    R[2,1] = z*y*(1-c) + x*s
    R[2,2] = c + (z**2)*(1-c)

    ret = np.dot(R,v)

    return ret

def matmul(A,v): #Matrix multiplication

    return np.array([np.dot(A[0],v),np.dot(A[1],v),np.dot(A[2],v)])

def beam_coor(beam,tele,ntele):
    # Converts beam coordinates from position to procetion on telescope
    #...check if tele is indeed a unit vector
    beam = -beam 
    zbeam = np.dot(beam,unit(tele))*unit(tele)
    ybeam = outplane(beam-zbeam,ntele)#...adjusted with unit 13-12-2018
    xbeam = beam-zbeam-ybeam
  
    return np.array([np.linalg.norm(xbeam),np.linalg.norm(ybeam),np.linalg.norm(zbeam)])
 
def tele_coor(tele_pos,beam,nbeam):
    # Converts tele coordinates from position beam components
    #...check if tele is indeed a unit vector
    ztele = np.dot(tele_pos,unit(beam))*unit(beam)
    ytele = outplane(tele_pos - ztele,nbeam)#...adjusted with unit 13-12-2018
    xtele = tele_pos-ztele-ytele
  
    return np.array([np.linalg.norm(xtele),np.linalg.norm(ytele),np.linalg.norm(ztele)])


def beam_ang(beam,tele,ntele):
    beam_conv = beam_coor(beam,unit(tele),ntele)
    mag = np.linalg.norm(beam_conv)
    ang_x = np.sin(abs(beam_conv[0])/mag)
    ang_y = np.sin(abs(beam_conv[1])/mag)

    return [ang_x,ang_y]
