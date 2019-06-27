from imports import *
import output
import numpy as np
import datetime
import os

def string_length(l,string):
    while len(string)<l:
        string = '0'+string

    return string

def get_date(option='date'):
    now = datetime.datetime.now()
    if option=='date':
        ret=string_length(2,str(now.year))+string_length(2,str(now.month))+string_length(2,str(now.day))
    elif option=='time':
        ret=string_length(2,str(now.hour))+string_length(2,str(now.minute))+string_length(2,str(now.second))
    #date=date+'-'+dir_extr
    return ret

def get_folder(direct=False,opt_date=True):
    if direct==False:
        if opt_date==True:
           date = get_date(option='date')+'/'
        elif opt_data==False:
            date==''
        direct = os.getcwd()+'/Results/'+date

    if not os.path.exists(direct):
        os.makedirs(direct)

    return direct

def write(inp,aim,title='',direct='',extr='',opt_date=True,opt_time=True,time='',extra_title='',include='all',exclude=[],offset=False,overwrite=True):
    
    date = get_date(option='date')
    if time=='':
        time = get_date(option='time')

    if direct=='':
        direct=get_folder(opt_date=opt_date)
    direct=direct+extr+'/'
    if not os.path.exists(direct):
        os.makedirs(direct)

    if opt_time==True:
        title=extra_title+'_'+time+'_'+title+'.txt'
    elif opt_time==False:
        title=extra_title+'_'+title+'.txt'

    name=direct+'/'+title
    check=name.split('.')
    if len(check)>=2:
        new_name=''
        for c in check[0:len(check)-1]:
            new_name=new_name+c+'.'
        name=new_name[0:len(new_name)-1]

    if overwrite==True and os.path.exists(name)==True:
        os.remove(name)

    writefile = open(name,'w')

    if offset!=False:
        writefile.write(str(offset))
        writefile.close()

        return direct+'/'+title
    
    else:
        writefile.write("BEGIN OPTIONS"+'\n')
        try:
            settings_all = [aim.data.stat,aim.aimset]
            for j in settings_all:
                for setting in j.__dict__.keys():
                    val = j.__dict__[setting]
                    write_on=True
                    try:
                        if '<' in str(val_new[v]):
                            write_on=False
                    except:
                        pass
                    if write_on:
                        writefile.write(setting+':: '+str(val)+':: '+str(type(val))+'\n')
        except AttributeError:
            pass
        writefile.write("END OPTIONS"+'\n')
        writefile.write('\n')

        for side in inp.__dict__.keys():
            for i in getattr(inp,side).__dict__.keys():
                for m in getattr(getattr(inp,side),i).__dict__.keys():
                    writefile.write('BEGIN\n')
                    outp = getattr(getattr(getattr(inp,side),i),m)
                    #writefile.write("Value:: "+str(m)+'\n')
                    writefile.write("SC:: "+str(i[-1])+'\n')
                    writefile.write("Side:: "+str(side)+'\n')
                    options = outp[-1]
                    options = options.split(', ')
                    for opt in options:
                        opt_split = opt.split('=')
                        writefile.write(opt_split[0]+':: '+opt_split[-1]+'\n')
                    for j in range(0,len(outp[0])):
                        values = outp[0][j]
                        values_str = '['
                        for v in values:
                            values_str = values_str+str(v)+','
                        values_str = values_str[0:-1] +']'
                        writefile.write(values_str+'\n')
                    writefile.write('END\n')
                    writefile.write('\n')
                                
        writefile.close()

        print(title+' saved in:')
        print(direct)

    return direct

def make_matrix(A):
    B = A.split(']')
    out=[]
    for j in B:
        if j!='':
            C = j.replace('[','')
            C = C.replace(']','')
            C = C.replace(',','')
            D = C.split(' ')
            try:
                D.remove('')
            except ValueError:
                pass
            clear=False
            while clear==False:
                try:
                    ind = D.index('')
                    del D[ind]
                except:
                    clear=True
                    pass
            E=[]
            #print(D)
            s = int(len(D)**0.5)
            for k in range(0,s):
                E.append(D[k*s:(k+1)*s])
            out.append(np.matrix(E,dtype=np.float64))
    
    return out

def read_options(filename,print_on=False,del_auto=True):
    aimset=utils.Object()
    read_on=False
    readfile = open(filename,'r')
    for line in readfile:
        if read_on:
            r=line
            rr=r.split(':: ')
            try:
                type_sel = rr[-1].split('\n')[0].split("'")[1]
                ret = getattr(np,type_sel)(rr[1])
                setattr(aimset,rr[0],ret)
            except:
                if type_sel=='numpy.float64':
                    setattr(aimset,rr[0],np.float64(rr[1]))
                if print_on:
                    print(rr[0],ret,type_sel)
                pass

        if 'BEGIN OPTIONS' in str(line):
            read_on=True
        if 'END OPTIONS' in str(line):
            break

    if del_auto==True:
        exceptions = ['dir_orbits','dir_savefig','filename','home','directory_imp','test_calc']
        for ex in exceptions:
            setattr(aimset, ex, getattr(settings.stat,ex))

    return aimset

def read_output(filenames=False,direct=False):
    ret=utils.Object()
    if type(direct)==str:
        direct=[direct]

    if filenames==False:
        f_list=[]
        for direct_calc in direct:
            for (dirpath, dirnames, filenames_calc) in os.walk(direct_calc):
                for f in filenames_calc:
                    if '.swp' not in f:
                        f_list.append(dirpath+'/'+f.split('/')[-1])

        filenames = f_list
    else:
        if type(filenames)!=list:
            filenames = [filenames]

    for filename in filenames:
        readfile = open(filename,'r')
        read_on=False
        for line in readfile:
            #print(line)
            if 'END\n'==line:
                if R['mode']=='mean_surface':

                    #B = R[R['value']][0][1]
                    R[R['value']][0][1] = make_matrix(R[R['value']][0][1])
                    pass
                R[R['value']][1] = 'value='+R['value']+', mode='+R['mode']
                try:
                    getattr(ret,R['Side'])
                except:
                    setattr(ret,R['Side'],utils.Object())
                    #getattr(Y,R['Side'])
                
                try:
                    getattr(getattr(ret,R['Side']),'i'+R['SC'])
                except:
                    setattr(getattr(ret,R['Side']),'i'+R['SC'],utils.Object())
                
                try: 
                    R_old = getattr(getattr(getattr(ret,R['Side']),'i'+R['SC']),R['value'])
                    #print(R['value'])
                    X = R[R['value']]
                    R_old[0][0] = R_old[0][0].append(X[0][0])
                    R_old[0][1] = R_old[0][1].append(X[0][1])
                    #print(R_old[1])
                    #print(X_[1])


                except AttributeError:
                    setattr(getattr(getattr(ret,R['Side']),'i'+R['SC']),R['value'],R[R['value']])
                
                #setattr(getattr(getattr(Y,R['Side']),'i'+R['SC']),[R['value']],R[R['value']])
                read_on=False
            elif read_on:
                if ':: ' in line:
                    [key, value] = line.split('\n')[0].split(':: ')
                    R[key] = value
                else:
                    if R['mode']=='mean_surface':
                        values=[]
                        value = line.split('\n')[0]
                        try:
                            R[R['value']]
                            new=False
                        except:
                            R[R['value']] = [[0,''],0]
                            new=True
                        if new==True:
                            R[R['value']][0][0] = np.array(utils.flatten(np.array(np.matrix(value))))
                        elif new==False:
                            R[R['value']][0][1] = R[R['value']][0][1]+value

                    else:
                        value = line.split('\n')[0]
                        value = value[1:-1]
                        values=[]
                        for v in value.split(','):
                            try:
                                values.append(np.float64(v))
                            except:
                                print(line)
                                print('Error')
                                print(R['mode'])

                        values = np.array(values)
                        try:
                            R[R['value']]
                        except:
                            R[R['value']] = [[],0]
                        R[R['value']][0].append(values)
                        
            elif 'BEGIN\n'==line:
                read_on=True
                R={}
        print(filename)
        readfile.close()
    options = read_options(filename)
    setattr(ret,'options',options)

    point_options = options.tele_control+'_'+options.PAAM_control+'__'+options.option_tele+'_'+options.option_PAAM

    point_settings = options.tele_method_solve+'_'+options.PAAM_method_solve+'__'+options.optimize_PAAM+'_'+str(options.optimize_PAAM_value).replace('.','d')

    try:
        getattr(ret,point_options)
    except:
        setattr(ret,point_options,utils.Object())

    setattr(getattr(ret,point_options),point_settings,ret)
    
    readfile.close()
    
    ret_new = ret
    del ret

    return ret_new,point_options,point_settings
