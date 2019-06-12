from pointLISA import *
import pointLISA
import os

rets = ['tele_ang','PAAM_ang','xoff','yoff','zoff','r','R','angx_ab_rec','angy_ab_rec','angx_nab_rec','angy_nab_rec','angx_ab_send','angy_ab_send','angx_nab_send','angy_nab_send','angx_rec','angy_rec','angx_send','any_send','angx_wf_send','angy_wf_send','piston','z_extra','FOV_wavefront','u','power']
#rets=[['xoff','yoff','zoff']]
SC=[1,2,3]
#SC=[1]
sides=['l','r']
#sides=['l']
#mode=['center','meanvar','mean_surface']
mode=['meanvar','mean_surface']
#mode=['center']
Nbins=2

def get_all(input_file=None,set_stat=utils.Object(),set_din=utils.Object()):
    data = pointLISA.run_stat.do_run(input_file=input_file,set_stat=set_stat)
    aim0= pointLISA.run_din.get_pointing(data,tele_control='no_control',PAAM_control='no_control',sampled=False,tele_ang_extra=False,PAAM_ang_extra=False,set_din=set_din,inp=False,init=True)
    aim= pointLISA.run_din.get_pointing(data,aim0=aim0,aim_old=aim0,sampled=False,set_din=set_din,init=False)

    return aim
folder='/home/ester/git/Results/20190606'
f_list=[]
for (dirpath, dirnames, filenames_calc) in os.walk(folder):
    for f in filenames_calc:
        if '.txt' in f:
            f_list.append(dirpath+'/'+f.split('/')[-1])

f_test=[f_list[0]]
for f in f_test:
    try:
        aim = get_all(input_file=f)
        
        try:
            t_plot
        except NameError:
            t0 = aim.data.t_all[3]
            tend = aim.data.t_all[-3]
            dt = aim.data.t_all[1] - aim.data.t_all[0]
            steps = 5
            t_plot = np.linspace(t0,tend,int((tend-t0)/dt)+1)

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
                        inp = out.make_functions(include=ret_inp,option='sampled',i=i,side=s,mode=m,t_plot=t_plot,Nbins=Nbins)
                        print(ret)
                        title = f.split('/')[-1]
                        title=title.split('_SC')[0]
                        #direct = f.split(title)[0]
                        direct = '/home/ester/git/Results/test/plots'
                        title=title+'_Res_'+tit_ret+'_SC_'+str(i)+'_side_'+s+'_mode_'+m
                        print(inp[1])
                        print('')
                        read_write.write(inp[1],aim,title=title,direct=direct,opt_date=False,opt_time=False,time='',extra_title='',include='all',exclude=[],offset=False)
                        print(direct+title)
    except:
        print(f+' has no proper information')
        #print(e)




