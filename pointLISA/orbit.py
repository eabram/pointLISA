##!/usr/bin/env python   
from pointLISA import *
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
	
    def make_function(self,**kwargs):
        return self.circular_orbit()

    def add_tidal(self,P,scale=1.0):
        cvec = lambda i,t: LA.unit(P(i,t))*(np.linalg.norm(P(i,t))-AU)
        scale = 1.0
        tidal_res = lambda i,t: scale*cvec(i,t)

        return lambda i,t: P(i,t) + tidal_res(i,t)

    def circular_orbit(self,lack=0.0,tidal=False,**kwargs):
        # AU is the distance between the LISA COM and the SUN
        # L is the LISA armlength
        # eta is the inclination angle
        # lack is what the initial angle of the LISA COM is (angle = tan(y/x))
        # T_SUN is the period of the LISA constellation orbiting the Sun
        # T_CW is the period of one cartwheel of LISA

        COM_x = lambda t: AU*np.cos((t/T_SUN)*2*np.pi-lack)
        COM_y = lambda t: AU*np.sin((t/T_SUN)*2*np.pi-lack)
        COM = lambda t: np.array([COM_x(t),COM_y(t),0])

        R = L_arm/(3.0**0.5)
        theta = lambda i,t: ((t/T_CW)*2*np.pi)+(i-1)*(2.0/3.0)*np.pi
        tri_x = lambda i,t: R*np.cos(theta(i,t))
        tri_y = lambda i,t: R*np.sin(theta(i,t))

        r_vec = lambda t: COM(t)/np.linalg.norm(COM(t))
        z_vec = lambda t: np.array([0,0,1])
        x_vec = lambda t: np.cross(r_vec(t),z_vec(t))


        y_vec = lambda t: LA.rotate(r_vec(t),x_vec(t),eta)

        x = lambda i,t: (tri_x(i,t)*x_vec(t))
        y = lambda i,t: (tri_y(i,t)*y_vec(t))

        P = lambda i,t: x(i,t)+y(i,t)+COM(t)

        if tidal==False:
            return P

        else:
            return self.add_tidal(P,scale=tidal)
  
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
            scale=np.float(self.scale)
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
            p = utils.calculations_constellation().high_precision(p_first)
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

