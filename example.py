import pointLISA
from pointLISA import *
import os
import matplotlib.pyplot as plt


length_calc=400 #days (can also be set to 'all')
relativistic=True

# Obtain STAT object (data, without pointing)
data = pointLISA.run_stat.do_run(length_calc=length_calc,relativistic=relativistic)

# Obtain pointing (aim)
aim0 = pointLISA.run_din.get_pointing(data,tele_control='no_control',PAAM_control='no_control') # initializing
aim = pointLISA.run_din.get_pointing(data,tele_control='full_control',PAAM_control='full_control',aim0=aim0,aim_old=aim0,option_tele='center',option_PAAM='center',sampled=True) # aim object with tele_control, PAAM_control, option_tele and option_PAAM setting. When sampled==True, a aim.aim_sampled attribute will be created. This aim object is sampled and fitted (iterpolated) for less accurate but fast follow up calculations

# Obtain output (out)
out_slow = pointLISA.OUTPUT(aim) # High precision, slow
out_fast = pointLISA.OUTPUT(aim.aim_sampled) # Low precision, fast









