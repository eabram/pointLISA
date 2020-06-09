##!/usr/bin/env python   
from pointLISA import *
# This class creates an orbit object which containes the positions of the spacecrafts ans its corresponding syntheticLISA object


class ORBIT():
    def __init__(self,input_param,**kwargs):
        for k in input_param.keys():
            setattr(self,k,input_param[k])

        if self.filename is not None:
            self.import_file(read_max=self.length_calc)
        elif self.filename is None:
            self.lisa_obj = self.circular_orbit(tidal=self.tidal)
            self.get_output()
        LISA,putp = self.get_LISA_object()
        self.LISA = LISA
        self.putp = putp
	
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
        if self.scale=='Default':
            print('Getting scale by filename:')
            a = self.filename
            a1 = a.split('.')[0]
            a1 = a1.split('_')
            for k in range(0,len(a1)):
                if 'scale' == a1[k]:
                    self.scale = float(a1[k+1])
            print(self.scale)
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
        print('')

        par=['t','p1x','p1y','p1z','p2x','p2y','p2z','p3x','p3y','p3z']
        p=[[],[],[]]
        t=[]
        filename=open(self.filename,'r')
        line_num=1
        scale=np.float(self.scale)
        line_count=0
        for line in filename.readlines():
            if read_max!='all':
                if line_count==read_max:
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
        p = self.high_precision(p_first)
        if self.timeunit == 'days':
            self.t=(np.array(t) - np.array([t[0]]*len(t)))*day2sec #... in sec, fist point at t=0
        elif self.timeunit == 'years':
            self.t=(np.array(t) - np.array([t[0]]*len(t)))*years2sec #... in sec, fist point at t=0
        else: #Already in seconds
            self.t=np.array(t) - np.array([t[0]]*len(t)) #... in sec, fist point at t=0
        Dt=self.t[1]-self.t[0] # Assuming Dt is constant
        if 'synthlisa' in str(type(self.interpolation_method)): #...kan weg
            self.lisa_obj = self.interpolation_method
        else:
            self.lisa_obj=SampledLISA(p[0],p[1],p[2],Dt,self.t[0],2) #1==lin_int, -1==lin_extr, 2==Lagrange, -2==derivative interpolator, 0==nearest
        self.p=p
        self.Dt=Dt
        self.pos=[self.t,p]
        self.par=par
        lisa_obj=self.lisa_obj
        self.linecount=line_count
        filename.close()

    def get_output(self):
        self.p=False
        self.Dt=False
        self.pos=False
        self.par=False
        self.linecount=False
        self.t = np.linspace(0,year2sec*10.0,int((year2sec*10.0)/day2sec))
    
    def high_precision(self,p):
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

    def fit_pointLISA(self):
        '''A function to fit the imported positional data'''
        def fit_twosteps(x,y):
            def fit_sin(x,a,b,c,d): #...Buiten functie zetten
                return a*np.sin(b*x+d)+c

            #Guesses
            pb = (2.0*np.pi)/(3600*24*365.25)
            pc = np.mean(y)
            pa = (np.max(y) - np.min(y))/2.0
            pd=0.0
            p0 = scipy.array([pa,pb,pc,pd])
            popt, pcov = scipy.optimize.curve_fit(fit_sin, x, y,p0=p0)
            y0 = np.array([fit_sin(t,popt[0],popt[1],popt[2],popt[3]) for t in x])
            f0 = lambda t: fit_sin(t,popt[0],popt[1],popt[2],popt[3])
            
            #y0 = [0]*len(y0)
            #f0 = lambda t: 0

            yrest = y-y0
            step=10
            select=5
            fits=[]
            f1_fits_tot = []
            starts=[]
            jstart=0
            residuals=[]
            rcond=[]
            while jstart<=len(x)-select:
                jend=jstart+step
                X = x[jstart:jend]
                Y = yrest[jstart:jend]
                order = [5]
                loop=True
                while loop is True:
                    polyfit0 = np.polyfit(X,Y,order[-1],full=True)
                    try:
                        if polyfit0[1][0]>2500:
                            order.append(order[-1] +1)
                            loop = True
                        else:
                            loop=False
                    except IndexError:
                        order.append(order[-1] - 1)
                        loop=True
                    try:
                        if order[-1]==order[-3]:
                            loop=False
                    except IndexError:
                        pass
                residuals.append(polyfit0[1])
                rcond.append(polyfit0[4])
                A = np.poly1d(polyfit0[0])
                fits.append(A)
                f1_fit = np.array([A(t) for t in X[0:select]])
                f1_fits_tot.append(f1_fit)
                starts.append(x[jstart])
                jstart = jstart+select

            def get_function(t,starts,fits):
                for i in range(1,len(starts)):
                    if starts[i]>t:
                        loc=i-1
                        break
                try:
                    return fits[loc](t)
                except UnboundLocalError:
                    #print(t)
                    return np.nan

            f1 = lambda t: get_function(t,starts,fits) +f0(t)

            return f1,residuals,rcond
        x = self.t
        F_all=[]
        residuals_all=[]
        rcond_all=[]
        for i in range(0,len(self.p)):
            F=[]
            residuals=[]
            rcond=[]
            for j in range(0,len(self.p[i][0])):
                y = self.p[i][:,j]
                calc = fit_twosteps(x,y)
                F.append(calc[0])
                residuals.append(calc[1])
                rcond.append(calc[2])
            F_all.append(F)
            residuals_all.append(residuals)
            rcond_all.append(rcond)

        putp = lambda i,t: np.array([F_all[i-1][0](t),F_all[i-1][1](t),F_all[i-1][2](t)])
        self.residuals = residuals_all
        self.rcond = rcond_all
        return putp

    def get_LISA_object(self):
        '''Creates the attribute LISA with a positional attribute putp'''
        if self.filename==None:
            lisa = utils.Object()
            lisa.putp = self
            putp = self
        else:
            type_select = self.interpolation_method
            print('type_select',type_select)
            if type_select=='sampled':
                lisa = self.lisa_obj
                putp = lambda i,t: np.array(lisa.putp(i,t))
            elif type_select=='py':
                func_nominal_arm = lambda i,time: nominal_arm(self,i,time) #...klopt niet
                lisa = PyLISA(self.lisa_obj,func_nominal_arm)
                putp = lambda i,t: np.array(lisa.putp(i,t))
            elif type_select=='cache':
                func_nominal_arm = lambda i,time: nominal_arm(self,i,time) #...klopt niet
                lisa_py = PyLISA(self.lisa_obj,func_nominal_arm)
                lisa = CacheLISA(lisa_py)
                putp = lambda i,t: np.array(lisa.putp(i,t))

            elif type_select=='interp1d':
                t_all = self.t
                pos = []
                x_interp = []
                y_interp = []
                z_interp = []
                for i in range(1,4):
                    x = self.p[i-1][:,0]
                    y = self.p[i-1][:,1]
                    z = self.p[i-1][:,2]

                    x_interp.append(calc.interpolate(t_all,x))
                    y_interp.append(calc.interpolate(t_all,y))
                    z_interp.append(calc.interpolate(t_all,z))

                putp = lambda i,t: np.array([x_interp[i-1](t),y_interp[i-1](t),z_interp[i-1](t)])
                lisa = utils.Object()
                #lisa.putp = putp
            elif type_select=='pointLISA':
                putp = self.fit_pointLISA()
                lisa = utils.Object()
                #lisa.putp = putp

        return lisa,putp
