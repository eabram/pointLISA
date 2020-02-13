from imports import *
import numpy as np
 
# This file contains general mathematical methodsi (linear algebra)

def norm(v):
    '''np.linalg.norm(v) function but shorter in notation'''
    return np.linalg.norm(v)

def unit(v):
    '''Returns the unit vector of v'''
    try:
        if norm(v)==0:
            return v #...adjust
            raise ValueError
        else:
            return v/norm(v)
    except:
        print('unit(v)')
        print(v)
        raise ValueError

def angle(v1,v2,dot=False):
    '''Calculates the angle between vector v1 and v2'''
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
        return np.nan

def inplane(v,n):
    '''Calculates the inplane component of v (with n the normal unit vector of its outplane)'''
    inplane_calc = v - (np.dot(v,n)/(np.linalg.norm(n)**2))*n
    return inplane_calc

def outplane(v,n):
    '''Calculates the outplane component of v (with n the normal unit vector of its outplane)'''
    outplane_calc = (np.dot(v,n)/(np.linalg.norm(n)**2))*n
    return outplane_calc

def ang_out(v,n):
    '''The angle of v with its outplane'''
    sign = np.sign(np.dot(outplane(v,n),n))
    return sign*angle(inplane(v,n),v)

def ang_in(v,n,r):
    '''The angle of v with its inplane'''
    inplane_calc = inplane(v,n)
    ang_in_calc = angle(inplane_calc,r)
    return ang_in_calc

def print_component(v,v_in,v_out,v_arm):
    '''Prints the normalized components'''
    n = norm(v)
    n_in = norm(v_in)
    n_out = norm(v_out)
    n_arm = norm(v_arm)

    print(n_in/n)
    print((n_out**2+n_in**2+n_arm**2)/n**2)
    print('')

    return 0

def ang_in_out(v1,v2,n,r,give='all'):
    '''Returns the inplane and/or outplane angle between v1 and v2 (with the same n and r vector)'''
    n = unit(n)
    v1_out = (np.dot(v1,n)*n)/(norm(n)**2)
    v2_out = (np.dot(v2,n)*n)/(norm(n)**2)

    ang_out_1 = np.arcsin(norm(v1_out)/norm(v1))
    ang_out_1 = ang_out_1 * np.sign(np.dot(v1_out,n))
    ang_out_2 = np.arcsin(norm(v2_out)/norm(v2))
    ang_out_2 = ang_out_2 * np.sign(np.dot(v2_out,n))

    v1_in = v1 - v1_out
    v2_in = v2 - v2_out

    ang_in_1 = angle(v1_in,r)
    ang_in_2 = angle(v2_in,r)
    ang_in = ang_in_1 - ang_in_2
    ang_out = ang_out_1 - ang_out_2

    if give=='all':
        return [ang_in,ang_out]
    elif give=='in':
        return ang_in
    elif give=='out':
        return ang_out

def rotate(v,n,ang,mag=False):
    '''Rotates v around n with angle ang'''
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
    '''Matrix multiplication'''
    return np.array([np.dot(A[0],v),np.dot(A[1],v),np.dot(A[2],v)])
