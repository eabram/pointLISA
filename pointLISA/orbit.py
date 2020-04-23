##!/usr/bin/env python  
from imports import *
import utils
from synthlisa import *
import numpy as np
import os
from fractions import Fraction
import math
import datetime
from decimal import *
from scipy.interpolate import interp1d
year2sec=32536000
day2sec = year2sec/365.25
 
# This class creates an orbit object which containes the positions of the spacecrafts ans its corresponding syntheticLISA object

class ORBIT():
    def __init__(self,input_param,**kwargs):
        for k in input_param.keys():
            setattr(self,k,input_param[k])
        if self.directory_imp != False:
            self.directory_imp=input_param['home']+self.directory_imp
        input_param['directory_imp'] =  self.directory_imp
        if self.filename=='None':
            print('Please select filename')
        else:
            self.import_file(read_max=self.length_calc)
        
    
    def import_file(self,read_max='all'):
        '''Import the file with the spacecrafts coordinates and time stamps and creates functions and a synthLISA object or in creates functions from a inported synthLISA object'''
        if 'function' not in str(type(self.LISA_opt)):
            directory=self.directory_imp
            file_orb=self.filename
             
            par=['t','p1x','p1y','p1z','p2x','p2y','p2z','p3x','p3y','p3z']
            p=[[],[],[]]
            t=[]
            direc=self.directory_imp
            if directory==False:
                file=open(file_orb,'r')
            else:
                file=open(direc+file_orb,'r')
            line_num=1
            scale=self.scale
            line_count=0
            for line in file.readlines():
                if read_max!='all':
                    if line_count==read_max+20:
                        break
                if line[0]!= '#':
                    a=line.split(' ')
                    cleanup=False
                    while cleanup==False:
                        try:
                            a.remove('')
                        except ValueError:
                            cleanup=True

                    b=[]
                    read_check=True
                    try:
                        for j in a:
                            b.append(np.float64(j))
                        a=b
                        p[0].append(np.array([a[1]*scale,a[2]*scale,a[3]*scale]))
                        p[1].append(np.array([a[4]*scale,a[5]*scale,a[6]*scale]))
                        p[2].append(np.array([a[7]*scale,a[8]*scale,a[9]*scale]))
                        t.append(a[0]) # ... [s]
                    except ValueError,e:
                        print(e)
                        print(a)
                        print()
                        read_check=False
                        pass
                    if read_check==True:
                        line_count=line_count+1
            p_first = np.array([p[0],p[1],p[2]],np.float64)
            p = utils.high_precision(p_first)
            if self.timeunit == 'days':
                self.t=(np.array(t) - np.array([t[0]]*len(t)))*day2sec #... in sec, fist point at t=0
            elif self.timeunit == 'years':
                self.t=(np.array(t) - np.array([t[0]]*len(t)))*years2sec #... in sec, fist point at t=0
            else: #Already in seconds
                self.t=np.array(t) - np.array([t[0]]*len(t)) #... in sec, fist point at t=0
            Dt=self.t[1]-self.t[0] # Assuming Dt s constant
            if 'synthlisa' in str(type(self.LISA_opt)):
                self.lisa_obj = self.LISA_opt
            else:
                self.lisa_obj=SampledLISA(p[0],p[1],p[2],Dt,self.t[0],2)
            self.p=p
            self.Dt=Dt
            self.pos=[self.t,p]
            self.par=par
            lisa_obj=self.lisa_obj
            self.linecount=line_count
        else:
            self.p=False
            self.Dt=False
            self.pos=False
            self.par=False
            self.linecount=False
            self.t = np.linspace(0,year2sec*10.0,int((year2sec*10.0)/day2sec))
