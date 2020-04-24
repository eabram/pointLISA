import __builtin__ as builtin
import datetime
from decimal import *
import scipy.fftpack
from fractions import Fraction
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

warnings.filterwarnings('ignore', 'The iteration is not making good progress')
warnings.simplefilter('ignore', np.RankWarning)

import parameters
para = parameters.__dict__
for k in para:
    globals()[k] = para[k]
import settings
import utils
for k in utils.__dict__.keys():
    if 'instance' in str(utils.__dict__[k]):
        if 'instance' in str(type(utils.__dict__[k])):
            globals()[k] = getattr(utils,k)
#import methods
import orbit
import static
import run_stat
#import LA

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
