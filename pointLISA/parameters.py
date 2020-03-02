from imports import *
import numpy as np
import os

year2sec=np.float64(32536000.0)
day2sec=np.float64(year2sec/365.25)
c=np.float64(300000000.0)

#Parameters:
labda =np.float64(1064.5*(10**-9)) # m, blz. 114
eta_opt = np.float64(0.23) # ...look up
eta_pd = np.float64(0.68) # A/W
P_L = np.float64(2) # W ...verify with new LISA technical speifications
P_min = 1.0e-12
#D = np.float64(0.20) # Diameter [m]
MAGNIFICATION = np.float64(135) # Check value in Technote

nu_0 = c/labda

h = np.float64(1.0/(6.241506*(10**18)))

#Parameters for point PAAM
PAA_out_lim = np.float64(0.5*0.000001)
PAA_out_marge = np.float64(0.1*0.000001)

D = np.float64(0.300) # blz. 12 Optical quality criterion of a truncated laser beam for extremely long distance heterodyne interferometry
gamma_0 = np.float64(1.12) # blz. 12 
L_tele = np.float64(1) # meter telecope length
E0 = np.float64(1) #...Adjust
FOV = np.float64(8e-6) #Field of View telescope
#Calculations
w0_laser = D/(2*gamma_0) # blz. 12
k = (2*np.pi)/labda
I_min = P_min/(((D**2)/4.0)*np.pi)
home_run = os.getcwd()



