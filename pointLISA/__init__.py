import __builtin__ as builtin
import copy
import datetime
from decimal import *
import scipy.fftpack
from fractions import Fraction
import inspect
import math
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
import scipy.optimize
from sympy import *
from synthlisa import *
import warnings
import yaml

warnings.filterwarnings('ignore', 'The iteration is not making good progress')
warnings.simplefilter('ignore', np.RankWarning)

#import parameters
direct = os.path.dirname(os.path.realpath(__file__))
filename = 'parameters.txt'
parfile = open(direct+'/'+filename,'r')
parameters_all = {}
for line in parfile:
    if '=' in line:
        [name,value] = line.split('=')
        name = name.replace(' ','')
        value.split('#')[0].replace(' ','')
        globals()[name] = eval(value)
        parameters_all[name] = globals()[name]
globals()['parameters_all'] = parameters_all
parfile.close()

import settings
import utils
for k in utils.__dict__.keys():
    if 'instance' in str(utils.__dict__[k]):
        if 'instance' in str(type(utils.__dict__[k])):
            globals()[k] = getattr(utils,k)
import orbit
import static
import run_stat
import output
import AIM
import run_din
