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

def write(inp,aimset,title='',direct='',extr='',opt_date=True,opt_time=True,time='',extra_title='',include='all',exclude=[]):
    
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
    writefile = open(direct+'/'+title,'w')

    #if len(inp)==1:
    #    inp=[inp]
    writefile.write("BEGIN OPTIONS"+'\n')
    for setting in aimset.__dict__.keys():
        val = aimset.__dict__[setting]
        write_on=True
        try:
            if '<' in str(val):
                write_on=False
        except:
            pass
        if write_on:
            writefile.write(setting+':: '+str(val)+':: '+str(type(val))+'\n')
    writefile.write("END OPTIONS"+'\n')
    writefile.write('\n')

    for side in inp.__dict__.keys():
        for i in getattr(inp,side).__dict__.keys():
            for m in getattr(getattr(inp,side),i).__dict__.keys():
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
                    writefile.write(str(outp[0][j])+'\n')
                writefile.write('\n')
                            
    writefile.close()

    print(title+' saved in:')
    print(direct)

    return direct


def read_options(filename):
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
                pass

        if 'BEGIN OPTIONS' in str(line):
            read_on=True
        if 'END OPTIONS' in str(line):
            break

    return aimset

