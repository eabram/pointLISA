#from __init__ import * 
from pointLISA import *
# This runner file imports a LISA data file, obtains a STAT object whith an ORBIT attribute 
def do_run(input_file=None,set_stat=utils.Object(),**kwargs):
    para = parameters.__dict__
    data_set0 = utils.get_settings(settings_input=input_file,select='stat')
        
    input_param={}
    for key,value in kwargs.items():
        input_param[key] = value
    
    for k in set_stat.__dict__.keys():
        if k not in input_param:
            input_param[k] = getattr(set_stat,k)

    for k in data_set0.__dict__.keys():
        if k not in input_param.keys():
            input_param[k] = data_set0.__dict__[k]

    for keys in input_param.keys():
        globals()[keys] = input_param[keys]

    filename_list=[]
    filename_done=[]
    PAA_res={}
    count=0

    for (dirpath, dirnames, filenames) in os.walk(dir_orbits):
        for i in filenames:
            if i.split('.')[-1]=='txt':
                a = dirpath+'/'+i
                a = a.replace('//','/')
                filename_list.append(a)

    data_all={}

    for i in filename_list:
        filename_name = i.split('/')[-1]
        
        if select == 'all':
            if '/try/' in i:
                execute = False
            else:
                execute = True
        else:
            if select in i:
                execute = True
            else:
                execute = False

        if filename_name in filename_done:
            execute = False

        if execute == True:
            filename_save = i.split('/')[-1].split('_')[0]
            data=static.STAT(input_param,para,filename = i).PAA_func()
            data.input_file = input_file

            filename_done.append(filename_name)
            count=count+1
            
            if test_calc==False: 
                data_all[filename_save] = data
                data_all[str(count)] = data
            
            print('test_calc: '+str(test_calc))
    
    if len(data_all.keys())==2:
        return data_all['1']
    else:
        return data_all
