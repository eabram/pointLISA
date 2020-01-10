from pointLISA import * 
import pointLISA
import os
import shutil
import copy

### This file imports poiting angles (and other properties) from a file and obtains a AIM object form it

### Settings ###
# Select datafiles for importing
folder_sel = "/20191017/_01/SS_2PAAM__20_days/"
folder_0='/home/ester/git/Results/'
source_folder = folder_0+folder_sel
folder=folder_0+folder_sel+'/data/'

# Select which properties will be ontained
rets=['tele_ang','PAAM_ang','adjust','offset']
#rets=[['xoff','yoff'],['I','I0','Ival','FOVlim'],['alpha','angx_wf_send','angy_wf_send','angx_wf_rec','angy_wf_rec']]

# Select for which arms and what kind of value (mode)
SC=[1,2,3]
sides=['l','r']
mode=['center']
Nbins=2
#################################################


clear_go=False
overwrite=True

def get_all(input_file=None,set_stat=utils.Object(),set_din=utils.Object()):
    options = pointLISA.read_write.read_options(input_file)
    set_stat_copy = utils.Object()
    for k in set_stat.__dict__.keys():
        setattr(set_stat_copy,k,getattr(set_stat,k))
    for k in pointLISA.settings.stat.__dict__.keys():
        if k not in set_stat.__dict__.keys():
            try:
                setattr(set_stat_copy,k,getattr(options,k))
            except AttributeError:
                print('No option for: '+k)
                pass
    data = pointLISA.run_stat.do_run(input_file=input_file,set_stat=set_stat_copy)

    set_din_copy = utils.Object()
    for k in set_din.__dict__.keys():
        setattr(set_din_copy,k,getattr(set_din,k))

    for k in pointLISA.settings.aimset.__dict__.keys():
        if k not in set_din_copy.__dict__.keys():
            try:
                setattr(set_din_copy,k,getattr(options,k))
            except AttributeError:
                print('No option for: '+k)
                pass

    aim= pointLISA.run_din.get_pointing(data,sampled=False,set_din=set_din_copy)
    return aim, set_stat_copy, set_din_copy

folder_sel = "/20191017/_01/SS_2PAAM__20_days/"
folder_0='/home/ester/git/Results/'
source_folder = folder_0+folder_sel
folder=folder_0+folder_sel+'/data/'

for (dirpath, dirnames, filenames_calc) in os.walk(source_folder):
    for f in filenames_calc:
        name=dirpath+'/'+f.split('/')[-1]
        new_name = name.split(source_folder)[-1]
        try:
            shutil.copy2(source_folder+new_name,folder+new_name)
        except:
            Q=new_name.split('/')
            q=''
            for i in range(0,len(Q)-1):
                q=q+'/'+Q[i]
            print(q)
            os.makedirs(folder+q+'/')
            print('')
            shutil.copy2(source_folder+new_name,folder+new_name)
 
run=True
get_output=True
if run==True:
    f_list=[]
    for (dirpath, dirnames, filenames_calc) in os.walk(folder):
        for f in filenames_calc:
            print(f)
            ext=False
            if '.' not in f:
                ext=True
            if ('.txt' in f) or ext==True:
                if 'Res' not in f:
                    f_list.append(dirpath+'/'+f.split('/')[-1])
                if 'Res' in f and clear_go==True:
                    os.remove(dirpath+'/'+f)
                    

    f_list.sort()
    f_test=f_list
    for f in f_test:
        #try:
        aim,set_stat_run, set_din_run = get_all(input_file=f)
        
        try:
            t_plot
        except NameError:
            t0 = aim.data.t_all[3]
            tend = aim.data.t_all[-3]
            dt = 300#aim.data.t_all[1] - aim.data.t_all[0]
            steps = 5
            t_plot = np.linspace(t0,tend,int((tend-t0)/dt)+1)
        
        if get_output==True:
            out = output.OUTPUT(aim)
            out.get_coordinate_systems(speed=1)
            for ret in rets:
                for i in SC:
                    for s in sides:
                        for m in mode:
                            if type(ret)!=list:
                                ret_inp=[ret]
                                tit_ret=ret
                            else:
                                ret_inp = ret
                                tit_ret=''
                                for tr in ret_inp:
                                    tit_ret = tit_ret+tr+'_'
                            print(ret_inp)

                            title = f.split('/')[-1]
                            title=title.split('_SC')[0]
                            direct = f.split(title)[0]
                            title=title+'_Res_'+tit_ret+'_SC_'+str(i)+'_side_'+s+'_mode_'+m+'.txt'
                            if (os.path.exists(title)==False) or (os.path.exists(title) and overwrite==True):
                                inp = out.make_functions(include=ret_inp,option='sampled',i=i,side=s,mode=m,t=t_plot,Nbins=Nbins)
                                print(ret)
                                print(inp[1])
                                print('')
                                read_write.write(inp[1],aim,title=title,direct=direct,opt_date=False,opt_time=False,time='',extra_title='',include='all',exclude=[],offset=False)
                                print(direct+title)
