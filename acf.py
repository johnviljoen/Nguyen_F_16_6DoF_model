#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 16:09:02 2020

@author: johnviljoen
"""

''' This file contains the aerodynamic coefficient functions (acf's) that run upon the lookup tables
created in preprocessing.py '''

import numpy as np
from scipy.interpolate import interpn, interp1d

# In[functions]

############## 3D ###############

def Cx(alpha,beta,dh):
    path_3D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/3D/'
    C_values_table = np.load(path_3D + 'Cx_values_table.npy')
    C_points_table = np.load(path_3D + 'Cx_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha, dh])
    C_lookup = interpn(C_points_table, C_values_table[:,:,:,0], point)
    return C_lookup[0]

def Cz(alpha,beta,dh):
    path_3D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/3D/'
    C_values_table = np.load(path_3D + 'Cz_values_table.npy')
    C_points_table = np.load(path_3D + 'Cz_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha, dh])
    C_lookup = interpn(C_points_table, C_values_table[:,:,:,0], point)
    return C_lookup[0]

def Cl(alpha,beta,dh):
    path_3D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/3D/'
    C_values_table = np.load(path_3D + 'Cl_values_table.npy')
    C_points_table = np.load(path_3D + 'Cl_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha, dh])
    C_lookup = interpn(C_points_table, C_values_table[:,:,:,0], point)
    return C_lookup[0]

def Cm(alpha,beta,dh):
    path_3D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/3D/'
    C_values_table = np.load(path_3D + 'Cm_values_table.npy')
    C_points_table = np.load(path_3D + 'Cm_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha, dh])
    C_lookup = interpn(C_points_table, C_values_table[:,:,:,0], point)
    return C_lookup[0]

def Cn(alpha,beta,dh):
    path_3D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/3D/'
    C_values_table = np.load(path_3D + 'Cn_values_table.npy')
    C_points_table = np.load(path_3D + 'Cn_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha, dh])
    C_lookup = interpn(C_points_table, C_values_table[:,:,:,0], point)
    return C_lookup[0]

################## 2D ##################

# Cy

def Cy(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cy_values_table.npy')
    C_points_table = np.load(path_2D + 'Cy_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

# lefs

def Cx_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cx_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cx_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

def Cy_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cy_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cy_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

def Cz_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cz_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cz_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

def Cl_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cl_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cl_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

def Cm_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cm_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cm_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

def Cn_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cn_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cn_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

# da = 20 deg

def Cy_a20(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cy_a20_values_table.npy')
    C_points_table = np.load(path_2D + 'Cy_a20_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

def Cn_a20(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cn_a20_values_table.npy')
    C_points_table = np.load(path_2D + 'Cn_a20_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

def Cl_a20(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cl_a20_values_table.npy')
    C_points_table = np.load(path_2D + 'Cl_a20_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

# da = 20 deg, lef

def Cy_a20_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cy_a20_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cy_a20_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

def Cn_a20_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cn_a20_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cn_a20_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

def Cl_a20_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cl_a20_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cl_a20_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

# r 30 deg

def Cy_r30(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cy_r30_values_table.npy')
    C_points_table = np.load(path_2D + 'Cy_r30_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

def Cn_r30(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cn_r30_values_table.npy')
    C_points_table = np.load(path_2D + 'Cn_r30_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

def Cl_r30(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cl_r30_values_table.npy')
    C_points_table = np.load(path_2D + 'Cl_r30_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup[0]

################# 1D ################

# p derivatives

def Cyp(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'Cyp_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'Cyp_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def Cnp(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'Cnp_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'Cnp_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def Clp(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'Clp_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'Clp_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def delta_Cyp_lef(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'delta_Cyp_lef_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'delta_Cyp_lef_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def delta_Cnp_lef(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'delta_Cnp_lef_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'delta_Cnp_lef_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def delta_Clp_lef(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'delta_Clp_lef_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'delta_Clp_lef_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

# q derivatives

def Cxq(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'Cxq_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'Cxq_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def Czq(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'Czq_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'Czq_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def Cmq(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'Cmq_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'Cmq_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def delta_Cxq_lef(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'delta_Cxq_lef_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'delta_Cxq_lef_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def delta_Czq_lef(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'delta_Czq_lef_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'delta_Czq_lef_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def delta_Cmq_lef(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'delta_Cmq_lef_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'delta_Cmq_lef_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

# r derivatives

def Cyr(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'Cyr_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'Cyr_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def Cnr(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'Cnr_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'Cnr_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def Clr(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'Clr_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'Clr_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def delta_Cyr_lef(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'delta_Cyr_lef_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'delta_Cyr_lef_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def delta_Cnr_lef(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'delta_Cnr_lef_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'delta_Cnr_lef_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def delta_Clr_lef(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'delta_Clr_lef_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'delta_Clr_lef_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

# other

def delta_Cm(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'delta_Cm_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'delta_Cm_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def delta_Cn_beta(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'delta_Cn_beta_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'delta_Cn_beta_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]

def delta_Cl_beta(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'delta_Cl_beta_values_table.npy')
    C_values_table = np.reshape(C_values_table,[C_values_table.size])
    C_points_table = np.load(path_1D + 'delta_Cl_beta_points_table.npy', allow_pickle=True)
    point = alpha
    f = interp1d(C_points_table, C_values_table)
    C_lookup = f(point)
    return C_lookup[0]