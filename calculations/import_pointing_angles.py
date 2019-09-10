from pointLISA import *
import pointLISA
import os
import shutil
import copy

#rets = ['tele_ang','PAAM_ang','xoff','yoff','zoff','r','R','angx_ab_rec','angy_ab_rec','angx_nab_rec','angy_nab_rec','angx_ab_send','angy_ab_send','angx_nab_send','angy_nab_send','angx_rec','angy_rec','angx_send','any_send','angx_wf_send','angy_wf_send','piston','z_extra','FOV_wavefront','u','power']
#rets = [['xoff','yoff','r']]
#rets=[['xoff','yoff','zoff','r'],['angx_ab_rec','angy_ab_rec','angx_ab_send','angy_ab_send'],['I','I0'],'FOV_wavefront',['piston','z_extra'],['angx_wf_send','angy_wf_send']]
#rets=[['angx_wf_rec','angy_wf_rec']]
rets=[['offset']]
#rets=[['xoff','yoff'],['angx_wf_send','angy_wf_send']]
#rets=[['PAAM_ang','tele_ang']]

SC=[1,2,3]
#SC=[1]
sides=['l','r']
#sides=['l']
#mode=['center','meanvar','mean_surface']
#mode=['meanvar','mean_surface']
mode=['center']
Nbins=2
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

    #aim0= pointLISA.run_din.get_pointing(data,tele_control='no_control',PAAM_control='no_control',sampled=False,tele_ang_extra=False,PAAM_ang_extra=False,set_din=set_din_copy,inp=False,init=True)
    #del set_din_copy
    #set_din_copy = copy.copy(set_din)
    #aim= pointLISA.run_din.get_pointing(data,aim0=aim0,aim_old=aim0,sampled=False,set_din=set_din_copy,init=False)
    aim= pointLISA.run_din.get_pointing(data,sampled=False,set_din=set_din_copy)
    return aim, set_stat_copy, set_din_copy

#folder_sel = '20190627/_02/'
#folder_sel = '20190626/_01/'
#folder_sel='20190810/_01/full_control__400_days/'#read_offset_ab_rel/'
#folder_sel='Data_for_Graphs/SSFCread/SS_full_control__400_days/read_offset_ab_rel/_option_tele_center_option_PAAM_center/'
folder_sel='/20190910/'
folder_0='/home/ester/git/Results/'
source_folder = folder_0+folder_sel
folder=folder_0+'test/'+folder_sel

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
    #f_test=[f_list[0]]
    f_test=f_list
    for f in f_test:
        #try:
        aim,set_stat_run, set_din_run = get_all(input_file=f)
        
        try:
            t_plot
        except NameError:
            t0 = aim.data.t_all[3]
            tend = aim.data.t_all[-3]
            dt = aim.data.t_all[1] - aim.data.t_all[0]
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
                                inp = out.make_functions(include=ret_inp,option='sampled',i=i,side=s,mode=m,t_plot=t_plot,Nbins=Nbins)
                                print(ret)
                                print(inp[1])
                                print('')
                                read_write.write(inp[1],aim,title=title,direct=direct,opt_date=False,opt_time=False,time='',extra_title='',include='all',exclude=[],offset=False)
                                print(direct+title)
                                #del inp




