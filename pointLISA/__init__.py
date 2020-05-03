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

import parameters
direct = os.path.dirname(os.path.realpath(__file__))
filename = 'parameters.txt'
parfile = open(direct+'/'+filename,'r')
parameter_names = []
for line in parfile:
    if '#' not in line:
        a = line.split('\n')[0].split(' ')[0]
        if len(a)>=1:
            parameter_names.append(a)
parfile.close()

for k in parameter_names:
    globals()[k] = getattr(parameters,k)
globals()['parameter_names'] = parameter_names

import orbit
import settings
import utils
for k in utils.__dict__.keys():
    if 'instance' in str(utils.__dict__[k]):
        if 'instance' in str(type(utils.__dict__[k])):
            globals()[k] = getattr(utils,k)
            print(k)
import orbit
import static
import run_stat
import output
import AIM
import run_din

#import yaml

