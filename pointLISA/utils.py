from pointLISA import *
import pointLISA
# This class contains different helper functions

#######################################################################

class Object(object):
    '''Creates a (empty) class'''
    pass

#######################################################################

class linear_algebra():
    # This class contains general mathematical methods (linear algebra)

    #def norm(self,v):
    #    '''np.linalg.norm(v) function but shorter in notation'''
    #    return np.linalg.norm(v)

    def unit(self,v):
        '''Returns the unit vector of v'''
        try:
            if np.linalg.norm(v)==0:
                return v #...adjust
                raise ValueError
            else:
                return v/np.linalg.norm(v)
        except:
            raise ValueError

    def angle(self,v1,v2,dot=False):
        '''Calculates the angle between vector v1 and v2'''
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        if norm_v1!=0 and norm_v2!=0:
            if dot==False:
                sin = np.linalg.norm(np.cross(v1,v2)/(norm_v1*norm_v2))
                return np.arcsin(sin)
            elif dot == True:
                cos = np.dot(v1,v2)/(norm_v1*norm_v2)
                return np.sign(np.dot(v1,v2))*np.arccos(cos)
        else:
            return np.nan
    
    def inplane_outplane(self,v,n):
        outplane_calc = (np.dot(v,n)/(np.linalg.norm(n)**2))*n
        inplane_calc = v - outplane_calc

        return inplane_calc.outplane_calc

    #def inplane(self,v,n):
    #    '''Calculates the inplane component of v (with n the normal unit vector of its outplane)'''
    #    inplane_calc = v - (np.dot(v,n)/(np.linalg.norm(n)**2))*n
    #    return inplane_calc

    #def outplane(self,v,n):
    #    '''Calculates the outplane component of v (with n the normal unit vector of its outplane)'''
    #    outplane_calc = (np.dot(v,n)/(np.linalg.norm(n)**2))*n
    #    return outplane_calc

    #def ang_out(self,v,n):
    #    '''The angle of v with its outplane'''
    #    sign = np.sign(np.dot(self.outplane(v,n),n))
    #    return sign*self.angle(self.inplane(v,n),v)

    #def ang_in(self,v,n,r):
    #    '''The angle of v with its inplane'''
    #    inplane_calc = self.inplane(v,n)
    #    ang_in_calc = self.angle(inplane_calc,r)
    #    return ang_in_calc

    #def print_component(self,v,v_in,v_out,v_arm):
    #    '''Prints the normalized components'''
    #    n = self.norm(v)
    #    n_in = self.norm(v_in)
    #    n_out = self.norm(v_out)
    #    n_arm = self.norm(v_arm)

    #    print(n_in/n)
    #    print((n_out**2+n_in**2+n_arm**2)/n**2)
    #    print('')

    #    return 0

    def ang_in_out_tot(self,v1,v2,n,r,give='tot'):
        '''Returns the inplane and/or outplane angle between v1 and v2 (with the same n and r vector)''' 
        if give=='tot':
            ang_tot = self.angle(v1,v2)
            return ang_tot
        else:
            n = self.unit(n)
            v1_out = (np.dot(v1,n)*n)/(np.linalg.norm(n)**2)
            v2_out = (np.dot(v2,n)*n)/(np.linalg.norm(n)**2)
            if give=='inp':
                v1_in = v1 - v1_out
                v2_in = v2 - v2_out
                ang_in_1 = self.angle(v1_in,r)
                ang_in_2 = self.angle(v2_in,r)
                ang_in = ang_in_1 - ang_in_2
                return ang_in
            elif give=='out':
                ang_out_1 = np.arcsin(np.linalg.norm(v1_out)/np.linalg.norm(v1))
                ang_out_1 = ang_out_1 * np.sign(np.dot(v1_out,n))
                ang_out_2 = np.arcsin(np.linalg.norm(v2_out)/np.linalg.norm(v2))
                ang_out_2 = ang_out_2 * np.sign(np.dot(v2_out,n))
                ang_out = ang_out_1 - ang_out_2
                return ang_out

    def rotate(self,v,n,ang,mag=False):
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

    def matmul(self,A,v): #Matrix multiplication
        '''Matrix multiplication'''
        return np.array([np.dot(A[0],v),np.dot(A[1],v),np.dot(A[2],v)])

#######################################################################

class calculations_constellation():
    def i_slr(self,i,side='all'):
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

    def get_settings(self,settings_input=None,select='stat',kwargs={}):
        ret = Object()
        original = Object()
        for k in settings.__dict__.keys():
            if select+'_' in k:
                setattr(original,k[len(select)+1:],settings.__dict__[k])
        for k in original.__dict__.keys():
            setattr(ret,k,original.__dict__[k])

        if settings_input!=None:
            setfile = open(settings_input,'r')
            for line in setfile:
                A = line.split(' ')
                name = A[0]
                value = A[-1].split('\n')[0]
                if name in ret.__dict__.keys():
                    delattr(ret,name)
                    try:
                        if value!='all':
                            try:
                                setattr(ret,name,eval(value))
                            except:
                                raise NameError()
                        else:
                            raise NameError('All is not a function')
                    except NameError:
                        setattr(ret,name,str(value))

        for key,value in kwargs.items():
            if key in ret.__dict__.keys():
                delattr(ret,key)
                setattr(ret,key,value)
            else:
                print('Option '+key+' is not a valid input parameter')

        return ret

#######################################################################

class calculations():
    #This class contains some (specific) calulation methods (the more general ones can be found in utils.py)
    def get_nearest_smaller_value(self,lst,val):
        '''Returns the nearest smaller vallue of val in list lst'''
        lst.sort()
        if val<lst[0]:
             pos = np.nan #...check if this holds
        else:
            for i in range(1,len(lst)):
                if val<lst[i] and val>=lst[i-1]:
                    pos = i-1
                    break

        try:
            return pos
        except UnboundLocalError:
            return np.nan
            pass

    def get_tele_SS(self,aim,method,i,t,side,x=False,y=False):
        '''Returns the pointing angle at time t for a SS control'''
        if method==False:
            if type(y)==bool:
                if side=='l':
                    fc = aim.tele_ang_l_fc
                elif side=='r':
                    fc = aim.tele_ang_r_fc
            else:
                fc=y
            t_adjust = x
            pos_t = self.get_nearest_smaller_value(t_adjust,t)
            
            if pos_t!=np.nan:
                if type(y)==bool:
                    try:
                        return fc(i,t_adjust[pos_t])
                    except:
                        #print(pos_t)
                        return np.nan
                else:
                    try:
                        return fc[pos_t]
                    except IndexError:
                        print(pos_t)
                        return np.nan

    def string_length(self,l,string):
        '''Returns the length of a string'''
        while len(string)<l:
            string = '0'+string

        return string

    def get_date(self,option='date'):
        '''Returns the date'''
        now = datetime.datetime.now()
        if option=='date':
            ret=self.string_length(2,str(now.year))+self.string_length(2,str(now.month))+self.string_length(2,str(now.day))
        elif option=='time':
            ret=self.string_length(2,str(now.hour))+self.string_length(2,str(now.minute))+self.string_length(2,str(now.second))
        return ret

    def get_folder(self,direct=False,opt_date=True):
        '''Returns a folder (name)'''
        if direct==False:
            if opt_date==True:
               date = self.get_date(option='date')+'/'
            elif opt_data==False:
                date==''
            direct = os.getcwd()+'/Results/'+date

        if not os.path.exists(direct):
            os.makedirs(direct)

        return direct

    def savefig(self,f,title='',direct=True,newtime=False,extension='.png'):
        '''This function can plot and save a figure'''
        if newtime==True:
            time = self.get_date(option='time')
        else:
            try:
                time
            except NameError:
                time='000000'
                pass
        
        date = self.get_date(option='date')

        if direct==True:
            direct = self.get_folder()
        
        if not os.path.exists(direct):
            os.makedirs(direct)
        
        title=direct+'/'+time+'-'+title+extension
        f.savefig(title)
        print('Saved as '+title)

        return 0

    def flatten(self,y):
        '''Returns a flattened list'''
        ynew=[]
        check=True
        try:
            len(y)
        except TypeError:
            ynew = [y]
            check=False
            pass

        if check==True:
            for i in range(0,len(y)):
                try:
                    for j in range(0,len(y[i])):
                        ynew.append(y[i][j])
                except TypeError:
                    ynew.append(y[i])

        return ynew

    def rdln(self,line,typ='text'):
        '''Reads information of a line from an imported file'''
        if '[array(' in line:
            newline = line.split('array(')
            line = newline[-1].split(')')[0]+']'
            A = line[0:-1]
            A = A.replace('[','')
            A = A.replace(']','')
            A = A.replace(' ','')
            A = A.split(',')
            B=[]
            for i in A:
                B.append(np.float64(i))
            B = B
            return [B]
        else:
            ret = line[0:-1]
            if typ=='float':
                return np.float64(ret)
            else:
                return ret

    def read(self,filename='',direct='',meas='all'):
        '''Reads imported values'''
        if type(meas)==str:
            meas = [meas]
        ret={}
        if direct=='':
            direct = self.get_folder()

        if filename=='':
            f_get=[]
            f_list=[]
            for (dirpath, dirnames, filenames) in os.walk(direct):
                #filenames.sort()
                for f in filenames:
                    f_list.append(dirpath+'/'+f.split('/')[-1])

            filenames=f_list
        else:
            print('Please select filename or leave blank')

        try:
            filenames
            go =True
        except UnboundLocalError:
            print('Please select proper title and/or directory')
            go=False
            pass

        if go==True:
            for filename_select in filenames:
                print('Reading '+filename_select)

                readfile = open(filename_select,'r')

                for line in readfile:
                    if 'Title' in line:
                        key1 = self.dln(line.split(':: ')[-1])
                        keys = self.rdln(line).replace(':',',').split(',')
                        print(keys)
                        key0 = (keys[3]+' ')[1:-1]
                        key1 = (keys[5]+' ')[1:-1]
                        if key0 not in ret.keys():
                            ret[key0] = {}
                        if key1 not in ret[key0].keys():
                            ret[key0][key1]={}
                    elif 'Iteration' in line:
                        iteration = self.rdln(line.split(':: ')[-1])
                        if iteration not in ret[key0][key1].keys():
                            ret[key0][key1][iteration] = {}
                    elif 'Option' in line:
                        option = self.rdln(line.split(':: ')[-1])
                        if option not in ret[key0][key1][iteration].keys():
                            ret[key0][key1][iteration][option]={}
                    elif 'ax_title' in line:
                        key2 = self.rdln(line.split(':: ')[-1])
                        if key2 not in ret[key0][key1][iteration][option].keys():
                            ret[key0][key1][iteration][option][key2]={}
                    elif 'Measurement' in line:
                        key2 = self.rdln(line.split(':: ')[-1])
                        if (key2.split(' ')[0] in meas) or (meas[0]=='all') and ('object' not in key2):
                            go=True
                            if key2 not in ret[key0][key1][iteration][option].keys():
                                ret[key0][key1][iteration][option][key2]={}
                        else:
                            go=False
                     
                    elif 'Label' in line:
                        if go==True:
                            key3 = self.rdln(line.split(':: ')[-1])
                            if key3 not in ret[key0][key1][iteration][option][key2].keys():
                                ret[key0][key1][iteration][option][key2][key3]={}
                                ret[key0][key1][iteration][option][key2][key3]['x']=np.array([])
                                ret[key0][key1][iteration][option][key2][key3]['y']=np.array([])

                    else:
                        if go==True:
                            try:
                                del x,y 
                            except NameError:
                                pass
                            try:
                                if ';' in line:
                                    [x,y] = line.split(';')
                                else:
                                    x = line
                                    y='np.nan'
                                ret[key0][key1][iteration][option][key2][key3]['x'] = np.append(ret[key0][key1][iteration][option][key2][key3]['x'],rdln(x,typ='float'))
                                try:
                                    ret[key0][key1][iteration][option][key2][key3]['y'] =np.append(ret[key0][key1][iteration][option][key2][key3]['y'],rdln(y,typ='float'))
                                    value=True
                                except ValueError:
                                    value=False
                            except ValueError,e:
                                print(e)
                                print(line)
                            if value==False:
                                ynew_list = self.rdln(y)[1:-1].split(' ')
                                ynew_write=[]
                                for ynew in ynew_list:
                                    try:
                                        ynew_write.append(np.float64(ynew))
                                    except:
                                        pass
                                ret[key0][key1][iteration][option][key2][key3]['y'] = np.append(ret[key0][key1][iteration][option][key2][key3]['y'],np.array(ynew_write))
                
                readfile.close()

        return ret



    ### Pointing functions

    ### Telescope pointing

    def get_wavefront_parallel(self,data,aim,i,t,side,PAAM_ang,ret,mode='opposite',precision=0,ksi=[0,0],angles=False):
        '''Calculates how the telescopes should rotate for a 90 degree angle between the recieving waveront and the receiving telescope'''
        [i_self,i_left,i_right] = const.i_slr(i)
        if mode=='opposite':
            if side=='l':
                tdel = data.L_sl_func_tot(i_self,t)
                if data.stat.calc_method=='Waluschka':
                    tdel0=tdel
                elif data.stat.calc_method=='Abram':
                    tdel0=0
                if angles==False:
                    tele_ang = aim.tele_l_ang(i_self,t+tdel0)
                else:
                    tele_ang=angles
                coor_start = beam_coor_out(data,i_self,t,tele_ang,PAAM_ang,aim.aimset.offset_tele['l'])
                coor_end = aim.tele_r_coor(i_left,t+tdel)
                start=aim.tele_l_start(i_self,t+tdel0)
                end=aim.tele_r_start(i_left,t+tdel)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]

            elif side=='r':
                tdel=data.L_sr_func_tot(i_self,t)
                if data.stat.calc_method=='Waluschka':
                    tdel0=tdel
                elif data.stat.calc_method=='Abram':
                    tdel0=0
                if angles==False:
                    tele_ang = aim.tele_r_ang(i_self,t+tdel0)
                else:
                    tele_ang=angles
                coor_start =  beam_coor_out(data,i_self,t,tele_ang,PAAM_ang,aim.aimset.offset_tele['r'])
                coor_end = aim.tele_l_coor(i_right,t+tdel)
                start = aim.tele_r_start(i_self,t+tdel0)
                end=aim.tele_l_start(i_right,t+tdel)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]

            [zoff,yoff,xoff]=LA.matmul(coor_start,end-start)
            if precision==0:
                R = zoff # Not precise
            elif precision==1:
                try:
                   [piston,z_extra] = wfe.z_solve(xoff,yoff,zoff,ret='all')
                except:
                    [piston,z_extra] = [np.nan,np.nan]
                R = wfe.R(piston)

            R_vec = np.array([(R**2-xoff**2-yoff**2)**0.5,yoff,xoff])
            tele_vec = LA.matmul(coor_start,-coor_end[0])
            angx_R = np.sign(R_vec[2])*abs(np.arctan(R_vec[2]/R_vec[0]))
            angy_R = np.sign(R_vec[1])*abs(np.arctan(R_vec[1]/R_vec[0]))
            angx_tele = np.sign(tele_vec[2])*abs(np.arctan(tele_vec[2]/tele_vec[0]))
            angy_tele = np.sign(tele_vec[1])*abs(np.arctan(tele_vec[1]/tele_vec[0]))
            angx = (angx_tele-angx_R)
            angy = (angy_tele-angy_R)
     
        elif mode=='self':
            if side=='l':
                tdel = data.L_rl_func_tot(i_self,t)
                if data.stat.calc_method=='Waluschka':
                    tdel0=tdel
                elif data.stat.calc_method=='Abram':
                    tdel0=0
              
                if angles==False:
                    tele_ang = aim.tele_r_ang(i_left,t-tdel)
                    tele_ang_end = aim.tele_l_ang(i_self,t-tdel0)
                    PAAM_ang = aim.beam_r_ang(i_left,t-tdel)
                elif len(angles)>=2:
                    tele_ang_end = angles[0]
                    tele_ang = angles[2]
                    PAAM_ang = aim.beam_r_ang(i_left,t-tdel)
                coor_start = beam_coor_out(data,i_left,t-tdel,tele_ang,PAAM_ang,aim.aimset.offset_tele['r'])
                coor_end = coor_tele(data,i_self,t,tele_ang_end)
                start = LA.unit(coor_start[0])*data.param.L_tele+data.putp(i_left,t-tdel)
                end = LA.unit(coor_end[0])*data.param.L_tele+data.putp(i_self,t-tdel0)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]
     
            elif side=='r':
                tdel = data.L_rr_func_tot(i_self,t)
                if data.stat.calc_method=='Waluschka':
                    tdel0=tdel
                elif data.stat.calc_method=='Abram':
                    tdel0=0

                if angles==False:
                    tele_ang = aim.tele_l_ang(i_right,t-tdel)
                    tele_ang_end = aim.tele_r_ang(i_self,t-tdel0)
                    PAAM_ang = aim.beam_l_ang(i_right,t-tdel)
                elif len(angles)>=2:
                    tele_ang_end = angles[0]
                    tele_ang = angles[2]
                    PAAM_ang = aim.beam_l_ang(i_right,t-tdel)
                coor_start = beam_coor_out(data,i_right,t-tdel,tele_ang,PAAM_ang,aim.aimset.offset_tele['l'])
                coor_end = coor_tele(data,i_self,t,tele_ang_end)
                start = LA.unit(coor_start[0])*data.param.L_tele+data.putp(i_right,t-tdel)
                end = LA.unit(coor_end[0])*data.param.L_tele+data.putp(i_self,t-tdel0)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]

                    
            [zoff,yoff,xoff]=LA.matmul(coor_start,end-start)
            out=OUTPUT(aim)

            if precision==0:
                R = zoff # Not precise
            elif precision==1:
                try:
                   [piston,z_extra] = out.z_solve(xoff,yoff,zoff,ret='all')
                except:
                    [piston,z_extra] = [np.nan,np.nan]
                R = out.R(piston)

            R_vec = np.array([(R**2-xoff**2-yoff**2)**0.5,yoff,xoff])
            R_vec_origin = LA.matmul(np.linalg.inv(coor_start),R_vec)
            R_vec_tele_rec = LA.matmul(coor_end,-R_vec_origin)
            angx = np.arctan(abs(R_vec_tele_rec[2]/R_vec_tele_rec[0]))*np.sign(R_vec_tele_rec[2])
            angy = np.arctan(abs(R_vec_tele_rec[1]/R_vec_tele_rec[0]))*np.sign(R_vec_tele_rec[1])

        if ret=='angy':
            return angy
        elif ret=='angx':
            return angx
        elif ret=='tilt':
            return (angx**2+angy**2)**0.5
        elif ret=='xoff':
            return xoff
        elif ret=='yoff':
            return yoff
        elif ret=='r':
            return (xoff**2 +yoff**2)**0.5

        elif ret=='all':
            ret_val={}
            ret_val['start']=start
            ret_val['end']=end
            ret_val['zoff']=zoff
            ret_val['yoff']=yoff
            ret_val['xoff']=xoff
            ret_val['coor_start']=coor_start
            ret_val['coor_end']=coor_end
            ret_val['bd_original_frame'] = np.array(coor_start[0])
            ret_val['bd_receiving_frame'] = LA.matmul(coor_end,ret_val['bd_original_frame'])
            ret_val['angx_func_rec'] = angx
            ret_val['angy_func_rec'] = angy
            ret_val['R_vec_tele_rec']=R_vec_tele_rec
            if precision==1:
                ret_val['piston']=piston
                ret_val['z_extra'] = z_extra
            ret_val['R']=R
            ret_val["R_vec_beam_send"] = R_vec
            ret_val['R_vec_origin'] = R_vec_origin
            ret_val['r']=(xoff**2+yoff**2)**0.5

            FOV_beamline = np.arccos(-ret_val['bd_receiving_frame'][0]/np.linalg.norm(ret_val['bd_receiving_frame']))
            FOV_wavefront = LA.angle(-R_vec_origin,coor_end[0])
            FOV_position = LA.angle(start-end,coor_end[0])
            ret_val['tilt']=FOV_wavefront
            ret_val['FOV_beamline']=FOV_beamline
            ret_val['FOV_wavefront']=FOV_wavefront
            ret_val['FOV_position']=FOV_position

            return ret_val

    def rotate_PAA_wavefront(self,data,aim,SC,t,side,ret,output_full=False):
        '''Rotates the telescope angles for a straignt hit wit the receiving wavefront'''
        [i_left,i_right,link] = const.i_slr(SC)

        f = lambda PAAM_ang,m: get_wavefront_parallel(data,aim,SC,t,side,PAAM_ang,m,mode='opposite',precision=0,ksi=[0,0],angles=False)
        ang_solve = scipy.optimize.brentq(lambda PAAM_ang: f(PAAM_ang,ret),np.float64(-0.1),np.float64(0.1))


        if output_full==True:
            return ang_solve,f(ang_solve,'yoff'),f(ang_solve,'angy')
        elif output_full==False:
            return ang_solve


    # Changes of coordinate system
    def coor_SC(self,data,i,t):
        '''Returns the coordinates of a spacecraft in [r,n,x] components'''
        t_calc=t

        r = LA.unit(data.r_func(i,t_calc))
        n = LA.unit(data.n_func(i,t_calc))
        x = np.cross(n,r)
        #offset = wfe.data.putp(i,t)

        return np.array([r,n,x])

    def coor_tele(self,data,i,t,ang_tele):
        '''Returns the coordinate system of telescope (same as SC but rotated over ang_tele inplane)'''
        L_tele = data.param.L_tele
        [r,n,x] = self.coor_SC(data,i,t)
        tele = r*L_tele
        tele = LA.rotate(tele,n,ang_tele)
        r = LA.unit(tele)
        x = np.cross(n,r)

        return np.array([r,n,x])

    def aberration_beam_coor(self,data,i,t,v,reverse=False): # if reverse==True: SUN-->SC, if reverse==False: SC-->SUN
        if data.stat.aberration==False:
            ret = v
        elif data.stat.aberration==True:
            V = data.vel.abs(i,t)
            if reverse==True:
                V=-V
            v_mag = np.linalg.norm(v)
            c_vec = LA.unit(v)*data.param.c
            ret = LA.unit(c_vec+V)*v_mag

        return ret

    def beam_coor_out__send(self,data,i,t,ang_tele,ang_paam,ang_tele_offset): # beam coordinates as seen from send frame, Sun coordinate
        '''Retunrs the coordinate system of the transmitted beam (same as SC but rotated over ang_tele inplane and ang_tele outplane)'''
        [r,n,x] = self.coor_tele(data,i,t,ang_tele+ang_tele_offset) #Telescope coordinate system

        r = LA.unit(LA.rotate(r,x,ang_paam)) # Rotate r in out of plane over ang_paam
        #r_new = aberration_beam_coor(data,i,t,r)
        n = np.cross(r,x)

        return np.array([r,n,x])

    def get_matrix_from_function(self,A,t):
        '''Returns a matrix from a function'''
        ret=[]
        for i in range(0,len(A)):
            vec=[]
            for j in range(0,len(A[i])):
                vec.append(A[i][j](t))
            ret.append(np.array(vec))

        return np.array(ret)

    def interpolate(self,x,y,method='interp1d'):
        '''Obtains a function by interpolation'''
        if method=='interp1d':
            if str(type(y[0]))!="<type 'numpy.ndarray'>":
                return interp1d(x,y,bounds_error=False)
            else:
                type_dim = str(type(y[0,0]))
                if type_dim!="<type 'numpy.ndarray'>":
                    ret=[]
                    for l in range(0,len(y[0])):
                        ret.append(interp1d(x,y[:,l],bounds_error=False))

                    return lambda t: np.array([ret[0](t),ret[1](t),ret[2](t)])
                else:
                    ret=[]
                    for i in range(0,len(y[0])):
                        vec=[]
                        for j in range(0,len(y[0][i])):
                            vec.append(interp1d(x,y[:,i,j],bounds_error=False))
                        ret.append(np.array(vec))
                    return lambda t: self.get_matrix_from_function(np.array(ret),t)
     
        else:
            print('Please select proper interpolation method (interp1d)')

    def SS_value(self,aim,link,t0,t_end,method,lim,ret='',tele_l=False,tele_r=False,option=False,print_on=False,value=0,offset_l=False,offset_r=False,dt=3600*100,scale=1): #set scale at maximum of <2.0
        '''Calculate the repointing time stamps and corresponfing telecsope angles'''

        if option==False:
            option = aim.aimset.option_tele
        
        if t_end>=aim.data.t_all[-1]:
            t_end = aim.data.t_all[-1]-dt

        tele_adjust_l = [] 
        tele_adjust_r = [] 
        offset_adjust_l = []
        offset_adjust_r = []
        out=OUTPUT(aim=aim)

        i = (link-2)%3

        [i_left,i_right,link] = const.i_slr(i)
     
        if ret=='Ivalx':
            lim = aim.data.param.P_min/(((aim.data.param.D**2)/4.0)*(np.pi))
            out_show = 'Ival'
        elif ret == 'angx_wf_rec':
            lim = aim.aimset.FOV
            out_show = 'alpha'

        print(lim)
        #t0 = aim.data.t_all[6]
        #t_end = t0+5*3600
        xtol = 1.0
        rtol = 1.0e-9
        step0 = 3600
        dt = 60.0
        print(t0,t_end)
        offset_l = False
        offset_r = False
        t_adjust = [t0]

        if aim.aimset.PAAM_deg==2:
            Done=False
            tele_l = aim.twoPAAM_tele_aim(i_left,t_adjust[-1],'l',test=True)[0]
            tele_r = aim.twoPAAM_tele_aim(i_right,t_adjust[-1],'r',test=True)[0]
            tele_adjust_l.append(tele_l)
            tele_adjust_r.append(tele_r)

            skip_l=False
            skip_r=False
            while Done==False:
                ang_l = lambda t: abs(aim.twoPAAM_tele_aim_SS_calc(i_left,t,'l',tele_adjust_l[-1])) - lim
                ang_r = lambda t: abs(aim.twoPAAM_tele_aim_SS_calc(i_right,t,'r',tele_adjust_r[-1])) - lim

                step=step0
                check=False
                while check==False and Done==False:
                    try:
                        t_l_new = scipy.optimize.brentq(ang_l,t_adjust[-1]+dt,t_adjust[-1]+step,xtol=xtol,rtol=rtol)
                        print(ang_l(t_l_new))
                        check=True
                        if t_adjust[-1]>t_end:
                            Done=True
                    except ValueError,e:
                        if str(e) =='f(a) and f(b) must have different signs':
                            if t_adjust[-1]+step>t_end:
                                Done=True
                                skip_l=True
                            else:
                                step = step*2
                                skip_l=False
                            pass

                step=step0
                check=False
                while check==False and Done==False:
                    try:
                        t_r_new = scipy.optimize.brentq(ang_r,t_adjust[-1]+dt,t_adjust[-1]+step,xtol=xtol,rtol=rtol)
                        print(ang_r(t_r_new))
                        check=True
                        if t_adjust[-1]>t_end:
                            Done=True
                    except ValueError,e:
                        if str(e) =='f(a) and f(b) must have different signs':
                            if t_adjust[-1]+step>t_end:
                                Done=True
                                skip_r=True
                            else:
                                step = step*2
                                skip_r=False
                            pass
                
                if Done==False:
                    write=True
                    if skip_l==True and skip_r==False:
                        t_adjust.append(t_r_new)
                    if skip_l==False and skip_r==True:
                        t_adjust.append(t_l_new)
                    if skip_l==False and skip_r==False:
                        t_adjust.append(np.minimum(t_l_new,t_r_new))
                    else:
                        write=False

                    if write==True:
                        tele_l = aim.twoPAAM_tele_aim(i_left,t_adjust[-1],'l',test=True)[0]
                        tele_r = aim.twoPAAM_tele_aim(i_right,t_adjust[-1],'r',test=True)[0]
                        tele_adjust_l.append(tele_l)
                        tele_adjust_r.append(tele_r)

                        print(t_adjust[-1]/t_end)

        if aim.aimset.PAAM_deg==1:
            Done=False

            if aim.aimset.option_tele=='wavefront':
                ret1 = 'Ival'
                print(ret,lim)
                [tele_ang_l_fc,tele_ang_r_fc] = aim.tele_control_ang_fc(option='wavefront',value=0.0)

                tele_l = tele_ang_l_fc(i_left,t_adjust[-1])
                tele_r = tele_ang_r_fc(i_right,t_adjust[-1])
                tele_adjust_l.append(tele_l)
                tele_adjust_r.append(tele_r)

                skip_l0=False
                skip_r0=False
                skip_l1=False
                skip_r1=False
                while Done==False:
                    send_l = lambda t: abs(getattr(self.values(aim,i_left,t,'l',tele_angle_l=tele_adjust_l[-1],tele_angle_r=tele_adjust_r[-1],beam_angle_l=False,beam_angle_r=False,offset_l=offset_l,offset_r=offset_r,ret=[ret],mode='rec'),ret)) - lim
                    send_lI = lambda t: getattr(self.values(aim,i_left,t,'l',tele_angle_l=tele_adjust_l[-1],tele_angle_r=tele_adjust_r[-1],beam_angle_l=False,beam_angle_r=False,offset_l=offset_l,offset_r=offset_r,ret=[ret1],mode='rec'),ret1) - aim.data.param.I_min

                    send_r = lambda t: abs(getattr(self.values(aim,i_right,t,'r',tele_angle_l=tele_adjust_l[-1],tele_angle_r=tele_adjust_r[-1],beam_angle_l=False,beam_angle_r=False,offset_l=offset_l,offset_r=offset_r,ret=[ret],mode='rec'),ret)) - lim
                    send_rI = lambda t: abs(getattr(self.values(aim,i_right,t,'r',tele_angle_l=tele_adjust_l[-1],tele_angle_r=tele_adjust_r[-1],beam_angle_l=False,beam_angle_r=False,offset_l=offset_l,offset_r=offset_r,ret=[ret1],mode='rec'),ret1)) - aim.data.param.I_min

                    step=step0
                    check=False
                    while check==False and Done==False:
                        try:
                            t_l_new0 = scipy.optimize.brentq(send_l,t_adjust[-1]+dt,t_adjust[-1]+step,xtol=xtol,rtol=rtol)
                            check=True
                            if t_adjust[-1]>t_end:
                                Done=True
                        except ValueError,e:
                            if str(e) =='f(a) and f(b) must have different signs':
                                if t_adjust[-1]+step>t_end:
                                    Done=True
                                    skip_l0=True
                                else:
                                    step = step*2
                                    skip_l0=False
                                pass

                    step=step0
                    check=False
                    while check==False and Done==False:
                        try:
                            t_l_new1 = scipy.optimize.brentq(send_lI,t_adjust[-1]+dt,t_adjust[-1]+step,xtol=xtol,rtol=rtol)
                            check=True
                            if t_adjust[-1]>t_end:
                                Done=True
                        except ValueError,e:
                            if str(e) =='f(a) and f(b) must have different signs':
                                if t_adjust[-1]+step>t_end:
                                    Done=True
                                    skip_l1=True
                                else:
                                    step = step*2
                                    skip_l1=False
                                pass
                    
                    t_l_new = min(t_l_new0,t_l_new1)
                    if skip_l0==True or skip_l1==True:
                        skip_l=True
                    else:
                        skip_l=False

                    step=step0
                    check=False
                    while check==False and Done==False:
                        try:
                            t_r_new0 = scipy.optimize.brentq(send_r,t_adjust[-1]+dt,t_adjust[-1]+step,xtol=xtol,rtol=rtol)
                            check=True
                            if t_adjust[-1]>t_end:
                                Done=True
                        except ValueError,e:
                            if str(e) =='f(a) and f(b) must have different signs':
                                if t_adjust[-1]+step>t_end:
                                    Done=True
                                    skip_r0=True
                                else:
                                    step = step*2
                                    skip_r0=False
                                pass
                    
                    step=step0
                    check=False
                    while check==False and Done==False:
                        try:
                            t_r_new1 = scipy.optimize.brentq(send_rI,t_adjust[-1]+dt,t_adjust[-1]+step,xtol=xtol,rtol=rtol)
                            check=True
                            if t_adjust[-1]>t_end:
                                Done=True
                        except ValueError,e:
                            if str(e) =='f(a) and f(b) must have different signs':
                                if t_adjust[-1]+step>t_end:
                                    Done=True
                                    skip_r1=True
                                else:
                                    step = step*2
                                    skip_r1=False
                                pass
                    
                    t_r_new = min(t_r_new0,t_r_new1)
                    if skip_r0==True or skip_r1==True:
                        skip_r=True
                    else: 
                        skip_r=False

                    if Done==False:
                        write=True
                        if skip_l==True and skip_r==False:
                            t_adjust.append(t_r_new)
                        if skip_l==False and skip_r==True:
                            t_adjust.append(t_l_new)
                        if skip_l==False and skip_r==False:
                            t_adjust.append(np.minimum(t_l_new,t_r_new))
                        else:
                            write=False
                        
                        if write==True:
                            A = self.tele_center_calc(aim,i_left,t_adjust[-1],scale=1,value=value,tele_l=None,tele_r=None,beam_l=None,beam_r=None,offset_l=False,offset_r=False)
                            [tele_l,tele_r] = A[0]
                            tele_adjust_l.append(tele_l)
                            tele_adjust_r.append(tele_r)

                            print(t_adjust[-1]/t_end)

            elif aim.aimset.option_tele == 'center':
                A = self.tele_center_calc(aim,i_left,t_adjust[-1],scale=1,value=value,tele_l=None,tele_r=None,beam_l=None,beam_r=None,offset_l=False,offset_r=False)
                [tele_l,tele_r] = A[0]
     
                #[tele_l,tele_r] = A[0]
                tele_adjust_l.append(tele_l)
                tele_adjust_r.append(tele_r)
                
                skip_l=False
                skip_r=False
                while Done==False:
                    send_l = lambda t: getattr(self.values(aim,i_left,t,'l',tele_angle_l=tele_adjust_l[-1],tele_angle_r=tele_adjust_r[-1],beam_angle_l=False,beam_angle_r=False,offset_l=offset_l,offset_r=offset_r,ret=[ret]),ret) - lim
                    send_r = lambda t: getattr(self.values(aim,i_right,t,'r',tele_angle_l=tele_adjust_l[-1],tele_angle_r=tele_adjust_r[-1],beam_angle_l=False,beam_angle_r=False,offset_l=offset_l,offset_r=offset_r,ret=[ret]),ret) - lim
                    
                    step=step0
                    check=False
                    while check==False and Done==False:
                        try:
                            t_l_new = scipy.optimize.brentq(send_l,t_adjust[-1]+dt,t_adjust[-1]+step,xtol=xtol,rtol=rtol)
                            check=True
                            if t_adjust[-1]>t_end:
                                Done=True
                        except ValueError,e:
                            if str(e) =='f(a) and f(b) must have different signs':
                                if t_adjust[-1]+step>t_end:
                                    Done=True
                                    skip_l=True
                                else:
                                    step = step*2
                                    skip_l=False
                                pass

                    step=step0
                    check=False
                    while check==False and Done==False:
                        try:
                            t_r_new = scipy.optimize.brentq(send_r,t_adjust[-1]+dt,t_adjust[-1]+step,xtol=xtol,rtol=rtol)
                            check=True
                            if t_adjust[-1]>t_end:
                                Done=True
                        except ValueError,e:
                            if str(e) =='f(a) and f(b) must have different signs':
                                if t_adjust[-1]+step>t_end:
                                    Done=True
                                    skip_r=True
                                else:
                                    step = step*2
                                    skip_r=False
                                pass

                    if Done==False:
                        write=True
                        if skip_l==True and skip_r==False:
                            t_adjust.append(t_r_new)
                        if skip_l==False and skip_r==True:
                            t_adjust.append(t_l_new)
                        if skip_l==False and skip_r==False:
                            t_adjust.append(np.minimum(t_l_new,t_r_new))
                        else:
                            write=False
                        
                        if write==True:
                            A = self.tele_center_calc(aim,i_left,t_adjust[-1],scale=1,value=value,tele_l=None,tele_r=None,beam_l=None,beam_r=None,offset_l=False,offset_r=False)
                            [tele_l,tele_r] = A[0]
                            tele_adjust_l.append(tele_l)
                            tele_adjust_r.append(tele_r)

                            print(t_adjust[-1]/t_end)
                
        return t_adjust,[tele_adjust_l,tele_adjust_r],i_left,i_right


    def tele_point_calc(self,aim,i,t,side,option,lim=False,method=False,value=0,scale=1,max_count=20,tele_l0=None,tele_r0=None,beam_l0=None,beam_r0=None,offset_l0=None,offset_r0=None,**kwargs): # Recommended to use aim0
        '''Calculates the (full control) telescope pointing angles (with the center or wavefront method)'''
        [i_self,i_left,i_right] = const.i_slr(i)
        if option=='center':
            if lim==False:
                lim = aim.aimset.limit_xoff
            if side=='l':
                ang = self.tele_center_calc(aim,i,t,lim=lim,value=value,tele_l=tele_l0,tele_r=tele_r0,beam_l=beam_l0,beam_r=beam_r0,offset_l=offset_l0,offset_r=offset_r0)[0][0]
            elif side=='r':
                ang = self.tele_center_calc(aim,const.i_slr(i)[2],t,lim=lim,value=value,tele_l=tele_l0,tele_r=tele_r0,beam_l=beam_l0,beam_r=beam_r0,offset_l=offset_l0,offset_r=offset_r0)[0][1]

        elif option=='wavefront':
            try:
                for k, value in kwargs.items:
                    locals()[k] = value
            except:
                pass
            if method==False:
                method = aim.aimset.tele_method_solve

            if lim==False:
                lim=aim.aimset.limit_angx

            if side=='l':
                ang = self.get_tele_wavefront(aim,i,t,'l',method,scale=scale,lim=lim,max_count=max_count,value=value,tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_l=beam_l0,beam_r=beam_r0,offset_l=offset_l0,offset_r=offset_r0)
            elif side=='r':
                ang = self.get_tele_wavefront(aim,i_right,t,'r',method,scale=scale,lim=lim,max_count=max_count,value=value,tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_l=beam_l0,beam_r=beam_r0,offset_l=offset_l0,offset_r=offset_r0)
                
        return ang


    def get_SS_func(self,x,y,x_check):
        '''Returns the SS function'''
        A = [t for t in x if t<x_check]
        val = y[len(A)-1]
        return np.float64(val)

    def t_sample(self,data,i,s,speed=1):
        '''Samples the timestamps'''
        if 'AIM' in str(data):
            data=data.data
        elif 'STAT' in str(data):
            pass
        else:
            raise(ValueError)
        
        t0 = data.t_all
        if speed==0:
            if s=='l':
                t_pref = np.array([t-data.L_rl_func_tot(i,t) for t in t0])
                t_next = np.array([t+data.L_sl_func_tot(i,t) for t in t0])
            elif s=='r':
                t_pref = np.array([t-data.L_rr_func_tot(i,t) for t in t0])
                t_next = np.array([t+data.L_sr_func_tot(i,t) for t in t0])
            
            t_sampled = np.concatenate((t0,t_pref))
            t_sampled = np.concatenate((t_sampled,t_next))
        elif speed==1:
            t_sampled = t0
        
        return np.sort(t_sampled)

    def get_t_sample(self,data,speed=0):
        '''Obtains the sampled timestamps'''
        t_l=[]
        t_r=[]
        t_all={}
        for i in range(1,4):
            t_l.append(t_sample(data,i,'l',speed=speed))
            t_r.append(t_sample(data,i,'r',speed=speed))
        
        t_all['l']= t_l
        t_all['r']= t_r

        return t_all




















    def get_coor_tele(self,aim,i,t,side,tele_angle=False):
        '''Gets telescope coordinate system'''
        if tele_angle==False:
            if side == 'l':
                try:
                    ret = aim.tele_l_coor(i,t)
                except AttributeError:
                    tele_angle = aim.tele_l_ang(i,t)
            elif side =='r':
                try:
                    ret = aim.tele_r_coor(i,t)
                except AttributeError:
                    tele_angle = aim.tele_r_ang(i,t)

        try:
            return ret
        except:
            ret = self.coor_tele(aim.data,i,t,tele_angle)
            return ret

    def get_coor_beam_in__sun(self,aim,i,t,tdel,side,tele_angle_send=False,beam_angle_send=False,tele_angle_rec=False,offset=False,out=3):
        '''Gets incoming (received) beam coordinate system'''
        [i_self,i_left,i_right] = const.i_slr(i)
        check=False
        if aim.data.stat.calc_method=='Abram':
            tdel0 = 0
        elif aim.data.stat.calc_method=='Waluschka':
            tdel0 = tdel
        try:
            if tele_angle_send==False and beam_angle_send==False:
                if side=='l':
                    u_sun = aim.beam_r_coor(i_left,t-tdel)
                elif side=='r':
                    u_sun = aim.beam_l_coor(i_right,t-tdel)
                check=True
        except AttributeError:
            check=False
            pass

        if check==False:
            if tele_angle_send is False:
                if side=='l':
                    tele_angle_send = np.radians(30.0)
                elif side=='r':
                    tele_angle_send = np.radians(-30.0)
            if beam_angle_send is False:
                if side=='l':
                    beam_angle_send = 0.0
                elif side=='r':
                    beam_angle_send = 0.0
            if offset is False:
                if side=='l':
                    offset = self.get_offset(aim,i_left,t-tdel,'r')
                elif side=='r':
                    offset = self.get_offset(aim,i_right,t-tdel,'l')
            elif offset == None:
                offset = 0.0

            if side=='l':
                u_sun = self.beam_coor_out(aim.data,i_left,t-tdel,tele_angle_send,beam_angle_send,offset)
            elif side=='r':
                u_sun = self.beam_coor_out(aim.data,i_right,t-tdel,tele_angle_send,beam_angle_send,offset)
        

        return u_sun

    def get_coor_beam_out__send(self,aim,i,t,side,tele_angle=False,beam_angle=False,offset=False):
        '''Gets outgoing (transmitted) beam coordinate system'''
        check=False
        
        if check==False:
            if tele_angle is False:
                if side == 'l':
                    tele_angle = aim.tele_l_ang(i,t)
                elif side =='r':
                    tele_angle = aim.tele_r_ang(i,t)
            elif tele_angle ==None:
                if side == 'l':
                    tele_angle = np.radians(-30.0)
                elif side =='r':
                    tele_angle = np.radians(30.0)

            if beam_angle is False:
                if side == 'l':
                    beam_angle = aim.beam_l_ang(i,t)
                elif side =='r':
                    beam_angle = aim.beam_r_ang(i,t)
            elif beam_angle==None:
                if side == 'l':
                    beam_angle = 0.0
                elif side =='r':
                    beam_angle = 0.0
            
            if offset is False:
                offset = self.get_offset(aim,i,t,side)
            elif offset==None:
                offset=0.0

            ret = self.beam_coor_out__send(aim.data,i,t,tele_angle,beam_angle,offset)
        return ret

    def get_offset(self,aim,i,t,side):
        '''Gets the offset angle (inplane) between the telescope and transmitted beam inplane angle'''
        try:
            ret = aim.offset[side][i](t)
        except TypeError:
            try:
                ret = aim.offset[side][i] 
            except TypeError:
                try:
                    ret = aim.offset(i,t,side)
                except TypeError:
                    print(i,t,side)
                    raise TypeError
        return ret

    def get_start_calc(self,aim,i,t,side,tele_angle):
        '''Gets the starting point of the telescope (where the tranmitted beam is leaving the telescope'''
        try:
            if side=='l':
                ret = aim.tele_l_start(i,t)
            elif side=='r':
                ret = aim.tele_r_start(i,t)
        except AttributeError:
            ret = np.array(aim.data.putp(i,t)) + LA.unit(self.get_coor_tele(aim,i,t,side,tele_angle=tele_angle)[0])*aim.data.param.L_tele

        return ret

    def values(self,inp,i,t,side,ksi=[0,0],mode='send',tele_angle_l=False,tele_angle_r=False,beam_angle_l=False,beam_angle_r=False,offset_l=False,offset_r=False,ret=[],getter=False):
        '''Runner function to obtain the output values for spacecraft i at time t'''
        [i_self,i_left,i_right] = const.i_slr(i)
        
        aim = inp
         
        if tele_angle_l==None:
            tele_angle_l = -np.radians(30.0)
        if tele_angle_r==None:
            tele_angle_r = np.radians(30.0)
        if beam_angle_l==None:
            beam_angle_l = 0.0
        if beam_angle_r==None:
            beam_angle_r = 0.0

        if getter is False:
            if mode=='send':
                if side=='l':
                    tdel = aim.data.L_sl_func_tot(i_self,t)
                    if aim.data.stat.calc_method=='Waluschka':
                        tdel0=tdel
                    elif aim.data.stat.calc_method=='Abram':
                        tdel0=0            
                    if offset_l is False:
                        offset_l = self.get_offset(aim,i_self,t+tdel0,'l')
                    elif offset_l == None:
                        offset_l = 0.0
                    if offset_r is False:
                        offset_r = self.get_offset(aim,i_left,t+tdel,'r')
                    elif offset_r == None:
                        offset_r = 0.0
                    i_send = i_self
                    i_rec = i_left

                elif side=='r':
                    tdel = aim.data.L_sr_func_tot(i_self,t)
                    if aim.data.stat.calc_method=='Waluschka':
                        tdel0=tdel
                    elif aim.data.stat.calc_method=='Abram':
                        tdel0=0
                    if offset_l is False:
                        offset_l = self.get_offset(aim,i_right,t+tdel,'l')
                    elif offset_l == None:
                        offset_l = 0.0
                    if offset_r is False:
                        offset_r = self.get_offset(aim,i_self,t+tdel0,'r')
                    elif offset_r == None:
                        offset_r = 0.0
                    i_send = i_self
                    i_rec = i_right

            elif mode=='rec':
                if side=='l':
                    tdel = aim.data.L_rl_func_tot(i_self,t)
                    if aim.data.stat.calc_method=='Waluschka':
                        tdel0=tdel
                    elif aim.data.stat.calc_method=='Abram':
                        tdel0=0
                    if offset_l is False:
                        offset_l = self.get_offset(aim,i_self,t-tdel0,'l')
                    elif offset_l==None:
                        offset_l=0.0
                    if offset_r is False:
                        offset_r = self.get_offset(aim,i_left,t-tdel,'r')
                    elif offset_r==None:
                        offset_r=0.0
                    i_send = i_left
                    i_rec = i_self

                elif side=='r':
                    tdel = aim.data.L_rr_func_tot(i_self,t)
                    if aim.data.stat.calc_method=='Waluschka':
                        tdel0=tdel
                    elif aim.data.stat.calc_method=='Abram':
                        tdel0=0
                    if offset_l is False:
                        offset_l = self.get_offset(aim,i_right,t-tdel,'l')
                    elif offset_l==None:
                        offset_l=0.0
                    if offset_r is False:
                        offset_r = self.get_offset(aim,i_self,t-tdel0,'r')
                    elif offset_r==None:
                        offset_r=0.0
                    i_send = i_right
                    i_rec = i_self

        elif getter is True:
            if mode=='send':
                i_send = i_self
                if side=='l':
                    i_rec = i_left
                    tdel = aim.data.L_sl_func_tot(i_self,t)
                    side_send = 'l'
                elif side=='r':
                    i_rec = i_right
                    tdel = aim.data.L_sr_func_tot(i_self,t)
                    side_send='r'
            elif mode=='rec':
                i_rec = i_self
                if side=='l':
                    i_send = i_left
                    tdel = aim.data.L_rl_func_tot(i_self,t)
                    side_send = 'r'
                elif side=='r':
                    i_send = i_right
                    side_send='l'
                    tdel = aim.data.L_rr_func_tot(i_self,t)
            if aim.data.stat.calc_method=='Waluschka':
                tdel0=tdel
            elif aim.data.stat.calc_method=='Abram':
                tdel0=0

        if getter is False:                
            if (mode=='send' and side=='l') or (mode=='rec' and side=='r'):
                tele_angle_start = tele_angle_l
                beam_angle_start = beam_angle_l
                tele_angle_end = tele_angle_r
                beam_angle_end = beam_angle_r
                offset_start = offset_l
                offset_end = offset_r
                
            elif (mode=='send' and side=='r') or (mode=='rec' and side=='l'):
                tele_angle_start = tele_angle_r
                beam_angle_start = beam_angle_r
                tele_angle_end = tele_angle_l
                beam_angle_end = beam_angle_l
                offset_start = offset_r
                offset_end = offset_l
    
        positions=Object()
        positions.method = aim.data.stat.calc_method
        positions.aim=aim
        positions.tele_angle_l = tele_angle_l
        positions.tele_angle_r = tele_angle_r
        positions.beam_angle_l = beam_angle_l
        positions.beam_angle_r = beam_angle_r
        positions.offset_l = offset_l
        positions.offset_r = offset_r
        
        param = ['mode','side','i_self','i_left','i_right','t','ksi','tdel','tdel0','aim','tele_angle_start','tele_angle_end','beam_angle_start','beam_angle_end','offset_l','offset_r','offset_start','offset_end','i_rec','i_send','getter']
        for p in param:
            try:
                setattr(positions,p,locals()[p])
            except Exception as e:
                #print(e)
                pass

        for r in ret:
            if r not in positions.__dict__.keys():
                #try:
                    #positions_new = getattr(outp,'get_'+r)(positions)
                positions_new = pointLISA.output.calc_value(positions,r)
                del positions
                positions = positions_new
                #except AttributeError,e:
                #    print(e)
                #    try:
                #        setattr(positions,r,getattr(aim,r)(i,t))
                #    except AttributeError,e:
                #        print(e)

        
        return positions

    def aberration(self,pos,u,mode='rec',**kwargs): # Only classical  
        if 'Object' in str(type(pos)):
            aber = pos.aim.data.stat.aberration
        elif 'STATIC' in str(pos):
            aber = pos.data.stat.aberration

        if aber==True:
            if 'Object' in str(type(pos)):
                if pos.aim.data.stat.aberration==False:
                    u_new = u
                elif pos.aim.data.stat.aberration==True:
                    if mode=='rec':
                        i = pos.i_rec
                        if pos.mode=='rec':
                            t = pos.t-pos.tdel0
                        elif pos.mode=='send':
                            t = pos.t-pos.tdel
                    elif mode=='send':
                        i = pos.i_send
                        if pos.mode=='rec':
                            t = pos.t+pos.tdel
                        elif pos.mode=='send':
                            t = pos.t+pos.tdel0
            elif 'STATIC' in str(pos):
                t = kwargs['t']
                i = kwargs['i']

            V = -pos.aim.data.vel.abs(i,t)
            u_mag = np.linalg.norm(u)
            c_vec = LA.unit(u)*c
            u_new = LA.unit(c_vec+V)*u_mag
        
        elif aber==False:
            u_new = u
        
        return u_new


    # Calculating the PAAM and telescope poiting angles
    def tele_center_calc(self,aim,i,t,scale=1,lim=1e-12,max_count=5,print_on=False,value=0,tele_l=False,tele_r=False,beam_l=False,beam_r=False,offset_l=False,offset_r=False):
        '''Obtains the telescope pointing angle when the telesope is pointed with the center method'''
        [i_self,i_left,i_right] = const.i_slr(i)
        
        lim = np.radians(5.0)
        if tele_l is False:
            tele_l=aim.tele_l_ang(i_self,t)
        elif tele_l==None:
            tele_l=np.radians(np.float64(-30.0))
        if tele_r is False:
            tele_r=aim.tele_r_ang(i_left,t)
        elif tele_r==None:
            tele_r=np.radians(np.float64(30.0))
        if beam_l is False:
            beam_l=aim.beam_l_ang(i_self,t)
        elif beam_l==None:
            beam_l=np.float64(0.0)
        if beam_r is False:
            beam_r=aim.beam_r_ang(i_self,t)
        elif beam_r==None:
            beam_r=np.float64(0.0)
        
        tele_l_old = tele_l
        tele_r_old = tele_r
        
        pos_send = lambda tele_l: self.values(aim,i_self,t,'l',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['xoff']).xoff
        send_solve = lambda tele_l: pos_send(tele_l)-value
        

        try:
            tele_l_new = scipy.optimize.brentq(send_solve,-lim-np.radians(30.0),lim-np.radians(30.0))
        except ValueError,e:
            if str(e)=='f(a) and f(b) must have different signs':
                tele_l_new=np.nan
     
        if tele_l_new!=np.nan:
            pos_rec = lambda tele_r: self.values(aim,i_left,t,'r',tele_angle_l=tele_l_new,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['xoff']).xoff
            rec_solve = lambda tele_r: pos_rec(tele_r)-value
            
            try:
                tele_r_new = scipy.optimize.brentq(rec_solve,-lim+np.radians(30.0),lim+np.radians(30.0))
            except ValueError,e:
                if str(e)=='f(a) and f(b) must have different signs':
                    tele_r_new=np.nan
        else:
            tele_r_new = np.nan
            
        return [[tele_l_new,tele_r_new], False]

    def PAAM_center_calc(self,aim,i,t,para='yoff_ab',scale=1,lim=1e-12,max_count=5,print_on=False,tele_l=None,tele_r=None,beam_l=None,beam_r=None,offset_l=None,offset_r=None,value=0,method='iter',margin=0.01):
        '''Obtains the PAAM pointing angle when the PAAM is pointed with the center method'''
        [i_self,i_left,i_right] = const.i_slr(i)
        
        if tele_l is False:
            tele_l = aim.tele_l_ang(i,t)
        elif tele_l==None:
            tele_l = np.radians(np.float64(-30.0))
        if tele_r is False:
            tele_r = aim.tele_r_ang(i,t)
        elif tele_r == None:
            tele_r = np.radians(np.float64(30.0))
        if beam_l is False:
            beam_l=aim.beam_l_ang(i,t)
        elif beam_l==None:
            beam_l = np.float64(0.0)
        if beam_r is False:
            beam_r=aim.beam_r_ang(i,t)
        elif beam_r==None:
            beam_r=np.float64(0.0)
        
        lim = 1.0e-3
        beam_l_old = beam_l
        beam_r_old = beam_r
        pos_send = lambda beam_l: self.values(aim,i_self,t,'l',tele_angle_l=tele_l,tele_angle_r=tele_l,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['yoff']).yoff
        send_solve = lambda beam_l: pos_send(beam_l)-value


        try:
            beam_l_new = scipy.optimize.brentq(send_solve,-lim,lim)
        except ValueError,e:
            if str(e)=='f(a) and f(b) must have different signs':
                beam_l_new=np.nan

        if beam_l_new!=np.nan:
            pos_rec = lambda beam_r: self.values(aim,i_left,t,'r',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam_l_new,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=['yoff']).yoff
            rec_solve = lambda beam_r: pos_rec(beam_r)-value

            try:
                beam_r_new = scipy.optimize.brentq(rec_solve,-lim,lim)
            except ValueError,e:
                if str(e)=='f(a) and f(b) must have different signs':
                    beam_r_new=np.nan
        else:
            beam_r_new = np.nan

        mode='Converged'

        return [[beam_l_new,beam_r_new], mode]

    def PAAM_wavefront_calc(self,aim,i,t,side,lim=1e-9,tele_l=False,tele_r=False):
        '''Obtains the PAAM pointing angle when the PAAM is pointed with the wavefront method'''
        if side=='l':
            angy = lambda beam: self.values(aim,i,t,'l',mode='send',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam,ret=['angy_wf_rec']).angy_wf_rec
        elif side=='r':
            angy = lambda beam: self.values(aim,i,t,'r',mode='send',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_r=beam,ret=['angy_wf_rec']).angy_wf_rec

        try:
            #ret = scipy.optimize.brentq(angy,-1e-5,1e-5,xtol=lim)
            ret = scipy.optimize.brentq(angy,-1e-1,1e-1,xtol=lim) #...solve between PAAM range (to do)
        except ValueError,e:
            if str(e)=='f(a) and f(b) must have different signs':
                ret=np.nan
        return ret



    def tele_wavefront_calc(self,aim,i_l,t,scale=1,lim=1e-12,max_count=5,print_on=False,value=0,tele_l=False,tele_r=False,beam_l=False,beam_r=False,offset_l=False,offset_r=False):
        '''Obtains the telescope pointing angle when the telesope is pointed with the center method'''
        [i_self,i_left,i_right] = const.i_slr(i_l)

        lim = np.radians(5.0)
        if tele_l is False:
            tele_l=aim.tele_l_ang(i_self,t)
        elif tele_l==None:
            tele_l=np.radians(np.float64(-30.0))
        if tele_r is False:
            tele_r=aim.tele_r_ang(i_left,t)
        elif tele_r==None:
            tele_r=np.radians(np.float64(30.0))
        if beam_l is False:
            beam_l=aim.beam_l_ang(i_self,t)
        elif beam_l==None:
            beam_l=np.float64(0.0)
        if beam_r is False:
            beam_r=aim.beam_r_ang(i_self,t)
        elif beam_r==None:
            beam_r=np.float64(0.0)
        
        tele_l_old = tele_l
        tele_r_old = tele_r
        
        par = 'angx_wf_rec'

        pos_l = getattr(self.values(aim,i_self,t,'l',mode='rec',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=[par]),par)
        tele_l_new = tele_l-pos_l
        
        pos_r = getattr(self.values(aim,i_left,t,'r',mode='rec',tele_angle_l=tele_l_new,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,offset_l=offset_l,offset_r=offset_r,ret=[par]),par)
        tele_r_new = tele_r-pos_r
        
        return [[tele_l_new,tele_r_new], [pos_l,pos_r]]


    def get_tele_wavefront(self,aim,i,t,side,method,scale=1,lim=1e-12,max_count=20,print_on=False,value=0.0,tele_angle_l=None,tele_angle_r=None,beam_l=None,beam_r=None,offset_l=False,offset_r=False): 
        '''Gets all telescope pointing angles along one arm for the wavefront method'''
        if side=='l':
            i_l = i
            tdel=0
        elif side=='r':
            i_l = const.i_slr(i)[2]
            i_r = i
            tdel = aim.data.L_rr_func_tot(i_r,t)
        
        tele_l_old = 0.0
        tele_r_old = 0.0
        if tele_angle_l==None:
            tele_l = np.radians(-30.0)
        else:
            tele_l = tele_angle_l
        if tele_angle_r==None:
            tele_r = np.radians(30.0)
        else:
            tele_r = tele_angle_r

        
        count=0
        while count<max_count:
            [[tele_l_new,tele_r_new],con] = self.tele_wavefront_calc(aim,i,t,tele_l=tele_l,tele_r=tele_r,beam_l=beam_l,beam_r=beam_r,offset_l=offset_l,offset_r=offset_r)
            count = count+1
            if count>= max_count:
                mode = 'Maximum iteration limit has been reached'
                tele_l = tele_l_new
                tele_r = tele_r_new           
                if print_on:
                    print(mode)
                break
            elif max(con)<1.0e-9: #max(tele_l_new-tele_l,tele_r_new-tele_r)<1.0e9:
                mode = 'Result is converged'
                tele_l = tele_l_new
                tele_r = tele_r_new
                if print_on:
                    print(mode)
                break
            tele_l = tele_l_new
            tele_r = tele_r_new
            
        if side=='l':
            return tele_l

        elif side=='r':
            return tele_r

#######################################################################

LA = linear_algebra()
const = calculations_constellation()
calc = calculations()
