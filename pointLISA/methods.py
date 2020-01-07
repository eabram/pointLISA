from imports import *   
import numpy as np
from scipy.interpolate import interp1d
from output import OUTPUT
import numpy as np
import output
import datetime
import os
from pointLISA import *
import scipy.optimize

#This file contains some (specific) calulation methods (the more general ones can be found in utils.py)

def get_putp_sampled(data,method='interp1d'):
    '''Returns an interpolation of the spacecraft positions'''
    t_all = data.orbit.t
    pos = []
    for i in range(1,4):
        pos_array=[]
        pos_x=[]
        pos_y=[]
        pos_z=[]
        for t in t_all:
            value = data.LISA.putp(i,t)
            if value[0]==0.0:
                pos_x.append(np.nan)
                pos_y.append(np.nan)
                pos_z.append(np.nan)
            else:
                pos_x.append(value[0])
                pos_y.append(value[1])
                pos_z.append(value[2])
        
        pos_x_interp  = interpolate(t_all,pos_x,method=method)
        pos_y_interp  = interpolate(t_all,pos_y,method=method)
        pos_z_interp  = interpolate(t_all,pos_z,method=method)
        pos.append([pos_x_interp,pos_y_interp,pos_z_interp])
        
    ret = lambda i,t: np.array([pos[i-1][0](t),pos[i-1][1](t),pos[i-1][2](t)])

    return ret


def get_nearest_smaller_value(lst,val):
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

def get_tele_SS(aim,method,i,t,side,x=False,y=False):
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
        pos_t = get_nearest_smaller_value(t_adjust,t)
        
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

def string_length(l,string):
    '''Returns the length of a string'''
    while len(string)<l:
        string = '0'+string

    return string

def get_date(option='date'):
    '''Returns the date'''
    now = datetime.datetime.now()
    if option=='date':
        ret=string_length(2,str(now.year))+string_length(2,str(now.month))+string_length(2,str(now.day))
    elif option=='time':
        ret=string_length(2,str(now.hour))+string_length(2,str(now.minute))+string_length(2,str(now.second))
    return ret

def get_folder(direct=False,opt_date=True):
    '''Returns a folder (name)'''
    if direct==False:
        if opt_date==True:
           date = get_date(option='date')+'/'
        elif opt_data==False:
            date==''
        direct = os.getcwd()+'/Results/'+date

    if not os.path.exists(direct):
        os.makedirs(direct)

    return direct

def savefig(f,title='',direct=True,newtime=False,extension='.png'):
    '''This function can plot and save a figure'''
    if newtime==True:
        time = get_date(option='time')
    else:
        try:
            time
        except NameError:
            time='000000'
            pass
    
    date = get_date(option='date')

    if direct==True:
        direct = get_folder()
    
    if not os.path.exists(direct):
        os.makedirs(direct)
    
    title=direct+'/'+time+'-'+title+extension
    f.savefig(title)
    print('Saved as '+title)

    return 0

def flatten(y):
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

def rdln(line,typ='text'):
    '''Reads information of a line from an imported file'''
    if '[array(' in line:
        newline = line.split('array(')
        line = newline[-1].split(')')[0]+']'
        A = line[0:-1]
        #print(A)
        #print('')
        A = A.replace('[','')
        A = A.replace(']','')
        A = A.replace(' ','')
        A = A.split(',')
        #print(A)
        B=[]
        for i in A:
            B.append(np.float64(i))
        B = B
        #print(B,len(B))
        return [B]
    else:
        ret = line[0:-1]
        if typ=='float':
            return np.float64(ret)
        else:
            return ret

def read(filename='',direct='',meas='all'):
    '''Reads imported values'''
    if type(meas)==str:
        meas = [meas]
    ret={}
    if direct=='':
        direct = get_folder()

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
            #print(filenames)
            print('Reading '+filename_select)

            readfile = open(filename_select,'r')

            for line in readfile:
                if 'Title' in line:
                    key1 = rdln(line.split(':: ')[-1])
                    keys = rdln(line).replace(':',',').split(',')
                    print(keys)
                    key0 = (keys[3]+' ')[1:-1]
                    key1 = (keys[5]+' ')[1:-1]
                    if key0 not in ret.keys():
                        ret[key0] = {}
                    if key1 not in ret[key0].keys():
                        ret[key0][key1]={}
                elif 'Iteration' in line:
                    iteration = rdln(line.split(':: ')[-1])
                    if iteration not in ret[key0][key1].keys():
                        ret[key0][key1][iteration] = {}
                elif 'Option' in line:
                    option = rdln(line.split(':: ')[-1])
                    if option not in ret[key0][key1][iteration].keys():
                        ret[key0][key1][iteration][option]={}
                elif 'ax_title' in line:
                    key2 = rdln(line.split(':: ')[-1])
                    if key2 not in ret[key0][key1][iteration][option].keys():
                        ret[key0][key1][iteration][option][key2]={}
                elif 'Measurement' in line:
                    key2 = rdln(line.split(':: ')[-1])
                    if (key2.split(' ')[0] in meas) or (meas[0]=='all') and ('object' not in key2):
                        go=True
                        if key2 not in ret[key0][key1][iteration][option].keys():
                            ret[key0][key1][iteration][option][key2]={}
                    else:
                        go=False
                 
                elif 'Label' in line:
                    if go==True:
                        key3 = rdln(line.split(':: ')[-1])
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
                            #except ValueError:
                            ynew_list = rdln(y)[1:-1].split(' ')
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

def get_wavefront_parallel(data,aim,i,t,side,PAAM_ang,ret,mode='opposite',precision=0,ksi=[0,0],angles=False):
    '''Calculates how the telescopes should rotate for a 90 degree angle between the recieving waveront and the receiving telescope'''
    [i_self,i_left,i_right] = utils.i_slr(i)
    if mode=='opposite':
        if side=='l':
            tdel = data.L_sl_func_tot(i_self,t)
            if data.calc_method=='Waluschka':
                tdel0=tdel
            elif data.calc_method=='Abram':
                tdel0=0
            if angles==False:
                tele_ang = aim.tele_l_ang(i_self,t+tdel0)
            else:
                tele_ang=angles
            coor_start = beam_coor_out(data,i_self,t,tele_ang,PAAM_ang,aim.offset_tele['l'])
            coor_end = aim.tele_r_coor(i_left,t+tdel)
            start=aim.tele_l_start(i_self,t+tdel0)
            end=aim.tele_r_start(i_left,t+tdel)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]

        elif side=='r':
            tdel=data.L_sr_func_tot(i_self,t)
            if data.calc_method=='Waluschka':
                tdel0=tdel
            elif data.calc_method=='Abram':
                tdel0=0
            if angles==False:
                tele_ang = aim.tele_r_ang(i_self,t+tdel0)
            else:
                tele_ang=angles
            coor_start =  beam_coor_out(data,i_self,t,tele_ang,PAAM_ang,aim.offset_tele['r'])
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
            if data.calc_method=='Waluschka':
                tdel0=tdel
            elif data.calc_method=='Abram':
                tdel0=0
          
            if angles==False:
                tele_ang = aim.tele_r_ang(i_left,t-tdel)
                tele_ang_end = aim.tele_l_ang(i_self,t-tdel0)
                PAAM_ang = aim.beam_r_ang(i_left,t-tdel)
            elif len(angles)>=2:
                tele_ang_end = angles[0]
                tele_ang = angles[2]
                PAAM_ang = aim.beam_r_ang(i_left,t-tdel)
            coor_start = beam_coor_out(data,i_left,t-tdel,tele_ang,PAAM_ang,aim.offset_tele['r'])
            coor_end = coor_tele(data,i_self,t,tele_ang_end)
            start = LA.unit(coor_start[0])*data.L_tele+data.putp(i_left,t-tdel)
            end = LA.unit(coor_end[0])*data.L_tele+data.putp(i_self,t-tdel0)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]

        
        elif side=='r':
            tdel = data.L_rr_func_tot(i_self,t)
            if data.calc_method=='Waluschka':
                tdel0=tdel
            elif data.calc_method=='Abram':
                tdel0=0

            if angles==False:
                tele_ang = aim.tele_l_ang(i_right,t-tdel)
                tele_ang_end = aim.tele_r_ang(i_self,t-tdel0)
                PAAM_ang = aim.beam_l_ang(i_right,t-tdel)
            elif len(angles)>=2:
                tele_ang_end = angles[0]
                tele_ang = angles[2]
                PAAM_ang = aim.beam_l_ang(i_right,t-tdel)
            coor_start = beam_coor_out(data,i_right,t-tdel,tele_ang,PAAM_ang,aim.offset_tele['l'])
            coor_end = coor_tele(data,i_self,t,tele_ang_end)
            start = LA.unit(coor_start[0])*data.L_tele+data.putp(i_right,t-tdel)
            end = LA.unit(coor_end[0])*data.L_tele+data.putp(i_self,t-tdel0)+coor_end[1]*ksi[1]+coor_end[2]*ksi[0]

                
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

def rotate_PAA_wavefront(data,aim,SC,t,side,ret,output_full=False):
    '''Rotates the telescope angles for a straignt hit wit the receiving wavefront'''
    [i_left,i_right,link] = utils.i_slr(SC)

    import scipy.optimize
    
    f = lambda PAAM_ang,m: get_wavefront_parallel(data,aim,SC,t,side,PAAM_ang,m,mode='opposite',precision=0,ksi=[0,0],angles=False)
    ang_solve = scipy.optimize.brentq(lambda PAAM_ang: f(PAAM_ang,ret),np.float64(-0.1),np.float64(0.1))


    if output_full==True:
        return ang_solve,f(ang_solve,'yoff'),f(ang_solve,'angy')
    elif output_full==False:
        return ang_solve


# Changes of coordinate system
def coor_SC(data,i,t):
    '''Returns the coordinates of a spacecraft in [r,n,x] components'''
    t_calc=t

    r = LA.unit(data.r_func(i,t_calc))
    n = LA.unit(data.n_func(i,t_calc))
    x = np.cross(n,r)
    #offset = wfe.data.putp(i,t)

    return np.array([r,n,x])

def coor_tele(data,i,t,ang_tele):
    '''Returns the coordinate system of telescope (same as SC but rotated over ang_tele inplane)'''
    L_tele = data.L_tele
    [r,n,x] = coor_SC(data,i,t)
    tele = r*L_tele
    tele = LA.rotate(tele,n,ang_tele)
    r = LA.unit(tele)
    x = np.cross(n,r)

    return np.array([r,n,x])

def beam_coor_out(data,i,t,ang_tele,ang_paam,ang_tele_offset):
    '''Retunrs the coordinate system of the transmitted beam (same as SC but rotated over ang_tele inplane and ang_tele outplane)'''
    [r,n,x] = coor_tele(data,i,t,ang_tele+ang_tele_offset) #Telescope coordinate system

    r = LA.unit(LA.rotate(r,x,ang_paam)) # Rotate r in out of plane over ang_paam
    n = np.cross(r,x)

    return np.array([r,n,x])

def i_slr(i):
    '''Returns [i_self,i_left,i_right]'''
    i_self = i
    i_left = (i+1)%3
    i_right = (i+2)%3

    i_ret = [i_self,i_left,i_right]
    for j in range(0,len(i_ret)):
        if i_ret[j]==0:
            i_ret[j]=3

    return i_ret

def get_matrix_from_function(A,t):
    '''Returns a matrix from a function'''
    ret=[]
    for i in range(0,len(A)):
        vec=[]
        for j in range(0,len(A[i])):
            vec.append(A[i][j](t))
        ret.append(np.array(vec))

    return np.array(ret)

def interpolate(x,y,method='interp1d'):
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
                return lambda t: get_matrix_from_function(np.array(ret),t)
 

    else:
        print('Please select proper interpolation method (interp1d)')

def SS_value(aim,link,t0,t_end,method,lim,ret='',tele_l=False,tele_r=False,option=False,print_on=False,value=0,offset_l=False,offset_r=False,dt=3600*24,scale=1): #set scale at maximum of <2.0
    '''Calculate the repointing time stamps and corresponfing telecsope angles'''
    import pointLISA

    if option==False:
        option = aim.aimset.option_tele
    
    if t_end>=aim.data.t_all[-1]:
        t_end = aim.data.t_all[-1]-dt

    t_adjust=[]
    tele_adjust_l = [] 
    tele_adjust_r = [] 
    offset_adjust_l = []
    offset_adjust_r = []

    i = (link-2)%3

    [i_left,i_right,link] = pointLISA.utils.i_slr(i)
    
    if method=='step':
        if tele_l ==False:
            tele_l = aim.tele_l_ang(i_left,t0)
        if tele_r==False:
            tele_r = aim.tele_r_ang(i_right,t0)
    
    t_adjust.append(t0)
    t_val  = t0

    if method=='step':
        while t_val<t_end:
            t_val = t_val+step

            func_l = getattr(output.values(aim,i_left,t_val,'l',ksi=[0,0],mode='send',tele_angle_l=tele_l,tele_angle_r=tele_r,ret=[ret]),ret)
            func_r = getattr(output.values(aim,i_right,t_val,'r',ksi=[0,0],mode='send',tele_angle_l=tele_l,tele_angle_r=tele_r,beam_angle_l=beam_l,beam_angle_r=beam_r,ret=[ret]),ret)

            if max(abs(func_l),abs(func_r))>=lim:
                t_adjust.append(t_val-step)
                tele_l = aim.tele_l_ang(i_left,t_val-step)
                tele_r = aim.tele_r_ang(i_right,t_val-step)
                offset_l = aim.offset['l'][i_left](t_val-step)
                offset_r = aim.offset['r'][i_right](t_val-step)
                tele_adjust_l.append(tele_l)
                tele_adjust_r.append(tele_r)
                offset_adjust_l.append(offset_l)
                offset_adjust_r.append(offset_r)

    elif method=='solve':
        out=output.OUTPUT(aim=aim)
        t_val = t_adjust[-1]

        if aim.PAAM_deg==1:
            tele_l = tele_point_calc(aim,i_left,t_val,'l',option,max_count=5,scale=1,value=value) 
            tele_r = tele_point_calc(aim,i_right,t_val,'r',option,max_count=5,scale=1,value=value) 
            offset_l = lambda t: False
            offset_r = lambda t: False

        elif aim.PAAM_deg==2:
            ret = 'ang_arm_tele_rec'
            ret_val = 'angx_arm_tele_rec'
            #A = aim.twoPAAM_pointing(i_left,t_val,'l',out,'rec')
            #B = aim.twoPAAM_pointing(i_right,t_val,'r',out,'rec')
            #tele_l = A[0]
            #tele_r = B[0]
            tele_l0 = np.radians(-30.0)
            tele_r0 = np.radians(30.0)
            offset_l=0.0
            offset_r=0.0
            beam_l=0.0
            beam_r=0.0
            f_l0 = lambda t, tele_l,tele_r: getattr(output.values(aim,i_left,t,'l',ksi=[0,0],mode='rec',tele_angle_l=tele_l,tele_angle_r=tele_r,offset_l=offset_l,offset_r=offset_r,beam_angle_l=beam_l,beam_angle_r=beam_r,ret=[ret]),ret)
            f_l1 = lambda t, tele_l,tele_r: getattr(output.values(aim,i_left,t,'l',ksi=[0,0],mode='rec',tele_angle_l=tele_l,tele_angle_r=tele_r,offset_l=offset_l,offset_r=offset_r,beam_angle_l=beam_l,beam_angle_r=beam_r,ret=[ret_val]),ret_val)
            f_r0 = lambda t, tele_l,tele_r: getattr(output.values(aim,i_right,t,'r',ksi=[0,0],mode='rec',tele_angle_l=tele_l,tele_angle_r=tele_r,offset_l=offset_l,offset_r=offset_r,beam_angle_l=beam_l,beam_angle_r=beam_r,ret=[ret]),ret)
            f_r1 = lambda t, tele_l,tele_r: getattr(output.values(aim,i_right,t,'r',ksi=[0,0],mode='rec',tele_angle_l=tele_l,tele_angle_r=tele_r,offset_l=offset_l,offset_r=offset_r,beam_angle_l=beam_l,beam_angle_r=beam_r,ret=[ret_val]),ret_val)
            tele_l_extra = f_l1(t_val,tele_l0,tele_r0)
            tele_r_extra = f_r1(t_val,tele_l0,tele_r0)
            tele_l = tele_l0+tele_l_extra
            tele_r = tele_r0+tele_r_extra
            tele_l_old = tele_l
            tele_r_old = tele_r

        tele_adjust_l.append(tele_l)
        tele_adjust_r.append(tele_r)
        
        end=False
        while (t_val<t_end) and end is False:
            if aim.PAAM_deg==1:
                f_l = lambda t: abs(getattr(output.values(aim,i_left,t,'l',ksi=[0,0],mode='rec',tele_angle_l=tele_l,tele_angle_r=tele_r,offset_l=offset_l(t),offset_r=offset_r(t),ret=[ret]),ret))-lim
                f_r = lambda t: abs(getattr(output.values(aim,i_right,t,'r',ksi=[0,0],mode='rec',tele_angle_l=tele_l,tele_angle_r=tele_r,offset_l=offset_l(t),offset_r=offset_r(t),ret=[ret]),ret))-lim
            elif aim.PAAM_deg==2:
                f_l = lambda t: abs(f_l0(t,tele_l,tele_r))-lim
                f_r = lambda t: abs(f_r0(t,tele_l,tele_r))-lim

            k=1
            found=False
            if t_val+dt*(k-1)>=t_end:
                #t_val=t_end
                end=True
                print('End of range')
                break

            while found==False and end is False:
                if t_val+dt*(k-1)>=t_end:
                    t_val=t_end
                    end=True
                    print('End of range')
                    break
                else:
                    under = t_val+dt*(k-1)+1.0
                    upper = t_val+dt*k
                    if under >t_end:
                        end=True
                        break
                    try:
                        t_l = scipy.optimize.brentq(f_l,under,upper,xtol=60.0)
                    except ValueError,e:
                        if str(e)=='f(a) and f(b) must have different signs':
                            t_l=np.inf
                            pass
                    try:
                        t_r = scipy.optimize.brentq(f_r,under,upper,xtol=60.0)
                    except ValueError,e:
                        if str(e)=='f(a) and f(b) must have different signs':
                            t_r=np.inf
                            pass
                    if t_l!=np.inf or t_r!=np.inf:
                        found=True
                    else:
                        k=k+1
         
            if found==True:
                t_adjust.append(min(t_l,t_r))
                t_val = t_adjust[-1]
                
                if aim.PAAM_deg==1:
                    tele_l = tele_point_calc(aim,i_left,t_val,'l',option,max_count=5,scale=1,value=value) 
                    tele_r = tele_point_calc(aim,i_right,t_val,'r',option,max_count=5,scale=1,value=value) 
                elif aim.PAAM_deg==2:
                    tele_l_extra = f_l1(t_val,tele_l0,tele_r0)
                    tele_r_extra = f_r1(t_val,tele_l0,tele_r0)
                    
                    sign_l = np.sign(tele_l_extra)
                    sign_r = np.sign(tele_r_extra)
                    #tele_l_new = tele_l0+tele_l_extra+(scale-1)*sign_l*aim.aimset.FOV
                    #tele_r_new = tele_r0+tele_r_extra+(scale-1)*sign_r*aim.aimset.FOV
                    
                    tele_l_new = tele_l0+tele_l_extra
                    tele_l_new = tele_l_new+(scale-1)*np.sign(tele_l_new-tele_l_old)*aim.aimset.FOV
                    tele_r_new = tele_r0+tele_r_extra
                    tele_r_new = tele_r_new+(scale-1)*np.sign(tele_r_new-tele_r_old)*aim.aimset.FOV


                    ##tele_l_new = tele_l_old+tele_l_extra*scale
                    ##tele_r_new = tele_r_old+tele_r_extra*scale
                    ##if scale!=1:
                    ##    c1 = f_l1(t_val,tele_l_new,tele_r_new)
                    ##    c2 = f_r1(t_val,tele_l_new,tele_r_new)
                    ##    if abs(c1)>lim:
                    ##        tele_l_new = tele_-(tele_l_new - tele_l_old) + tele_l_old
                    ##
                    ##    if abs(c2)>lim:
                    ##        tele_r = -(tele_r_new - tele_r_old) + tele_r_old


                    #
                    #tele_l_new = tele_l0+tele_l_extra
                    #tele_r_new = tele_r0+tele_r_extra

                    #tele_l = (tele_l_new - tele_l_old)*scale + tele_l_old
                    #tele_r = (tele_r_new - tele_r_old)*scale + tele_r_old
                    tele_l = tele_l_new
                    tele_r = tele_r_new
                    tele_l_old = tele_l_new
                    tele_r_old = tele_r_new


                tele_adjust_l.append(tele_l)
                tele_adjust_r.append(tele_r)

                if print_on==True:
                    print(t_val/t_end,t_val,tele_l_extra,tele_r_extra)
                    #print(t_end/(3600*24.0))
                    print(f_l1(t_val,tele_l,tele_r),f_r1(t_val,tele_l,tele_r))
                    print(f_l1(t_val,tele_adjust_l[-2],tele_adjust_r[-2]),f_r1(t_val,tele_adjust_l[-2],tele_adjust_r[-2]))

            
    return t_adjust,[tele_adjust_l,tele_adjust_r],i_left,i_right


def tele_point_calc(aim,i,t,side,option,lim=False,method=False,value=0,scale=1,max_count=20,tele_l0=None,tele_r0=None,beam_l0=None,beam_r0=None,offset_l0=None,offset_r0=None,**kwargs): # Recommended to use aim0
    '''Calculates the (full control) telescope pointing angles (with the center or wavefront method)'''
    if option=='center':
        if lim==False:
            lim = aim.aimset.limit_xoff
        if side=='l':
            ang = output.tele_center_calc(aim,i,t,lim=lim,value=value,tele_l=tele_l0,tele_r=tele_r0,beam_l=beam_l0,beam_r=beam_r0,offset_l=offset_l0,offset_r=offset_r0)[0][0]
        elif side=='r':
            ang = output.tele_center_calc(aim,utils.i_slr(i)[2],t,lim=lim,value=value,tele_l=tele_l0,tele_r=tele_r0,beam_l=beam_l0,beam_r=beam_r0,offset_l=offset_l0,offset_r=offset_r0)[0][1]

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
            ang = output.get_tele_wavefront(aim,i,t,'l',method,scale=scale,lim=lim,max_count=max_count,value=value,tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_l=beam_l0,beam_r=beam_r0,offset_l=offset_l0,offset_r=offset_r0)
        elif side=='r':
            ang = output.get_tele_wavefront(aim,i,t,'r',method,scale=scale,lim=lim,max_count=max_count,value=value,tele_angle_l=tele_l0,tele_angle_r=tele_r0,beam_l=beam_l0,beam_r=beam_r0,offset_l=offset_l0,offset_r=offset_r0)

    return ang


def get_SS_func(x,y,x_check):
    '''Returns the SS function'''
    A = [t for t in x if t<x_check]
    val = y[len(A)-1]
    return np.float64(val)

def t_sample(data,i,s,speed=1):
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

def get_t_sample(data,speed=0):
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
