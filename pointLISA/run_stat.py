from imports import *
import os
#import plotfile2
#import save_fig
#import writefile

def do_run(input_param={},**kwargs):
    

    para = parameters.__dict__
    setting = settings.__dict__
    for key,value in kwargs.items():
        input_param[key] = value
        print(key,value)
        #else:
        #    print(input_param)
        #    raise ValueError("Double input variable for "+key)

    for k in setting.keys():
        if '__'  in k:
            del setting[k]
        else:
            if k not in input_param.keys():
                input_param[k] = setting[k]
                #globals()[k] = settings[k]
    
    for keys in input_param.keys():
        globals()[keys] = input_param[keys]


    filename_list=[]
    filename_done=[]
    PAA_res={}
    count=0

    for (dirpath, dirnames, filenames) in os.walk(dir_orbits):
        print(filenames)
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
        print('Dir_extr:'+dir_extr)
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
            data=STAT(input_param,para,filename = i).PAA_func() 
            filename_done.append(filename_name)
            count=count+1

            if test_calc==False:
                if plot_on==True:
                    data = plotfile2.do_plot(data,dir_extr,i,new_folder,tstep,plot_on=plot_on)
                    data = writefile.do_writefile(data,data_use=True)
                    save_fig.do_save_fig(data)
                
                data_all[filename_save] = data
                data_all[str(count)] = data
            
            print('test_calc: '+str(test_calc))
    
    return data_all
