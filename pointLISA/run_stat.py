from imports import * 
import os
#import plotfile2
#import save_fig
#import writefile

def do_run(input_file=None,input_param={},set_stat=utils.Object,**kwargs):
    para = parameters.__dict__
    setting = settings.stat.__dict__
    if input_file!=None:
        setting_new = read_write.read_options(input_file).__dict__
        for k in setting.keys():
            try:
                setting[k] = setting_new[k]
            except KeyError:
                pass

    
    input_param_new={}
    for key,value in kwargs.items():
        input_param_new[key] = value
    
    for k in set_stat.__dict__.keys():
        if k not in input_param_new:
            input_param_new[k] = getattr(set_stat,k)

    for k in setting.keys():
        if '__'  in k:
            del setting[k]
        else:
            if k not in input_param_new.keys():
                input_param_new[k] = setting[k]

    del input_param
    input_param = input_param_new
    
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
        if i == filename_list[0]:
            new_folder=False # Adjust if you (don't) want to override
        else:
            new_folder=False
        
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
            #for k in input_param.keys():
            #    print(k,input_param[k])
            data=STAT(input_param,para,filename = i).PAA_func()
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
