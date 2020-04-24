#from imports import *

import numpy as np
import os
from fractions import Fraction
import math
import datetime
from decimal import *
from scipy.interpolate import interp1d

import parameters
import orbit
import static
import run_stat


para = parameters.__dict__
for k in para:
    globals()[k] = para[k]

#import os
#global os
#import time
#global time
#import yaml
#global yaml
#
#from synthlisa import *
#global np
#import numpy as np
#import matplotlib.pyplot as plt
#from fractions import Fraction
#import math
#import datetime
#from scipy.interpolate import interp1d
#from scipy.interpolate import RegularGridInterpolator
#import scipy.optimize
#from sympy import * 

#import orbit
#import read_write
#import run_stat
#import settings
#import static
#import utils





# Read parameters
#para = parameters.__dict__
#for k in para:
#    globals()[k] = para[k]
#print('Done initialization of package')
