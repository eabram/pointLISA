from pointLISA import *

### This file will read and/or write output values (from and to datafiles) and also contains some helper functions
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

def write(inp,title='',direct='',extr='',opt_date=True,opt_time=True,time='',extra_title='',include='all',exclude=[],offset=False,overwrite=True): #Done
    '''Writes the output to a datafile'''

    aim = inp.options['aim']
    
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
    
    if title[-4:len(title)]!='.txt':
        title=title+'.txt'



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
            settings_all = [aim.data.constellationset,aim.alignmentset]
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
        
        t_all = inp.options['t_all']
        for side in inp.__dict__.keys():
            if side in ['l','r']:
                for i in getattr(inp,side).__dict__.keys():
                    if i[0:2]=='SC':
                        print(i)
                        #for m in getattr(getattr(inp,side),i).keys():
                        for m in inp.options['cases']:
                            writefile.write('BEGIN\n')
                            outp = getattr(getattr(inp,side),i)[m]
                            writefile.write("SC:: "+str(i[-1])+'\n')
                            writefile.write("Side:: "+str(side)+'\n')
                            writefile.write("mode:: "+str(inp.options['mode'])+'\n')
                            writefile.write("case:: "+str(m)+'\n')
                            values_str = '['
                            for j in range(0,len(t_all)):
                                value = t_all[j]
                                if j==len(t_all)-1:
                                    values_str = values_str + str(value)+']'
                                else:
                                    values_str = values_str+str(value)+','
                            writefile.write(values_str+'\n')
                            values_str = '['
                            for j in range(0,len(outp)):
                                value = outp[j]
                                if j==len(outp)-1:
                                    values_str = values_str + str(value)+']'
                                else:
                                    values_str = values_str+str(value)+','
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

def read_options(filename,print_on=False,del_auto=True): #Done
    '''Reads the options (used settings) from the datafile filenme'''
    aimset=utils.Object()
    read_on=False
    readfile = open(filename,'r')
    for line in readfile:
        if read_on:
            r=line
            rr=r.split(':: ')
            try:
                type_sel = rr[-1].split('\n')[0].split("'")[1]
                if type_sel=='bool':
                    if 'False' in rr[1]:
                        ret = False
                    elif 'True' in rr[1]:
                        ret =True
                else:
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

    aimset_ex=utils.Object()
    if del_auto==True:
        exceptions = ['dir_orbits','dir_savefig','filename','home','directory_imp','test_calc']
        for ex in exceptions:
            try:
                setattr(aimset_ex, ex, getattr(settings.stat,ex))
            except AttributeError:
                pass

    aimset_ret = utils.Object()
    for k in aimset_ex.__dict__.keys():
        setattr(aimset_ret,k,getattr(aimset_ex,k))
    for k in aimset.__dict__.keys():
        if k not in aimset_ret.__dict__.keys():
            setattr(aimset_ret,k,getattr(aimset,k))

    return aimset_ret

def read_output(filenames=False,direct=False):
    '''Reads the properties (output) of filename'''
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
        ret_values=utils.Object()
        for line in readfile:
            #print(line)
            if 'END\n'==line:                
                print(R['value'])
                try:
                    getattr(ret_values,R['Side'])
                except AttributeError:
                    setattr(ret_values,R['Side'],utils.Object())
                try:
                    getattr(getattr(ret_values,R['Side']),'SC'+R['SC'])
                except AttributeError:
                    setattr(getattr(ret_values,R['Side']),'SC'+R['SC'],{})
                getattr(getattr(ret_values,R['Side']),'SC'+R['SC'])[R['case']] = R['value']

                read_on=False
            elif read_on:
                if ':: ' in line:
                    [key, value] = line.split('\n')[0].split(':: ')
                    R[key] = value
                    print(key,value)
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
                        R['value'].append(values)
                    except:
                        R['value'] = []
                        R['value'].append(values)
                    if len(R['value'])==2:
                        print(R['value'])

                        
            elif 'BEGIN\n'==line:
                read_on=True
                R={}
        print(filename)
        readfile.close()
    options = read_options(filename)
    setattr(ret_values,'options',options)
    #setattr(ret,'options',values)

    #point_options = options.tele_control+'_'+options.PAAM_control+'__'+options.option_tele+'_'+options.option_PAAM

    #point_settings = options.tele_method_solve+'_'+options.PAAM_method_solve+'__'+options.optimize_PAAM+'_'+str(options.optimize_PAAM_value).replace('.','d')

    #try:
    #    getattr(ret,point_options)
    #except:
    #    setattr(ret,point_options,utils.Object())

    #setattr(getattr(ret,point_options),point_settings,ret)
    
    readfile.close()
    

    return ret_values#,point_options,point_settings
