#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 16:09:02 2020

@author: johnviljoen
"""

''' This file contains the aerodynamic coefficient functions (acf's) that run upon the lookup tables
created in preprocessing.py '''

import numpy as np
from scipy.interpolate import interpn
from sys import exit


# In[functions]

alpha_in = 90
beta_in = -30
dh_in = -25

point = np.array([beta_in,alpha_in,dh_in])

############## 3D ###############

def Cx(alpha,beta,dh):
    path_3D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/3D/'
    C_values_table = np.load(path_3D + 'Cx_values_table.npy')
    C_points_table = np.load(path_3D + 'Cx_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha, dh])
    C_lookup = interpn(C_points_table, C_values_table[:,:,:,0], point)
    return C_lookup

def Cz(alpha,beta,dh):
    path_3D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/3D/'
    C_values_table = np.load(path_3D + 'Cz_values_table.npy')
    C_points_table = np.load(path_3D + 'Cz_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha, dh])
    C_lookup = interpn(C_points_table, C_values_table[:,:,:,0], point)
    return C_lookup

def Cl(alpha,beta,dh):
    path_3D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/3D/'
    C_values_table = np.load(path_3D + 'Cl_values_table.npy')
    C_points_table = np.load(path_3D + 'Cl_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha, dh])
    C_lookup = interpn(C_points_table, C_values_table[:,:,:,0], point)
    return C_lookup

def Cm(alpha,beta,dh):
    path_3D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/3D/'
    C_values_table = np.load(path_3D + 'Cm_values_table.npy')
    C_points_table = np.load(path_3D + 'Cm_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha, dh])
    C_lookup = interpn(C_points_table, C_values_table[:,:,:,0], point)
    return C_lookup

def Cn(alpha,beta,dh):
    path_3D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/3D/'
    C_values_table = np.load(path_3D + 'Cn_values_table.npy')
    C_points_table = np.load(path_3D + 'Cn_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha, dh])
    C_lookup = interpn(C_points_table, C_values_table[:,:,:,0], point)
    return C_lookup

################## 2D ##################

# Cy

def Cy(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cy_values_table.npy')
    C_points_table = np.load(path_2D + 'Cy_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

# lefs

def Cx_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cx_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cx_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

def Cy_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cy_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cy_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

def Cz_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cz_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cz_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

def Cl_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cl_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cl_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

def Cm_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cm_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cm_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

def Cn_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cn_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cn_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

# da = 20 deg

def Cy_a20(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cy_a20_values_table.npy')
    C_points_table = np.load(path_2D + 'Cy_a20_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

def Cn_a20(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cn_a20_values_table.npy')
    C_points_table = np.load(path_2D + 'Cn_a20_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

def Cl_a20(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cl_a20_values_table.npy')
    C_points_table = np.load(path_2D + 'Cl_a20_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

# da = 20 deg, lef

def Cy_a20_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cy_a20_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cy_a20_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

def Cn_a20_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cn_a20_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cn_a20_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

def Cl_a20_lef(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cl_a20_lef_values_table.npy')
    C_points_table = np.load(path_2D + 'Cl_a20_lef_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

# r 30 deg

def Cy_r30(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cy_r30_values_table.npy')
    C_points_table = np.load(path_2D + 'Cy_r30_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

def Cn_r30(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cn_r30_values_table.npy')
    C_points_table = np.load(path_2D + 'Cn_r30_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

def Cl_r30(alpha,beta):
    path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'
    C_values_table = np.load(path_2D + 'Cl_r30_values_table.npy')
    C_points_table = np.load(path_2D + 'Cl_r30_points_table.npy', allow_pickle=True)
    point = np.array([beta, alpha])
    C_lookup = interpn(C_points_table, C_values_table[:,:,0], point)
    return C_lookup

################# 1D ################

# p derivatives

def Cxq(alpha):
    path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'
    C_values_table = np.load(path_1D + 'Cxq_values_table.npy')
    C_values_table = np.reshape(C_values_table,[111])
    print(C_values_table.shape)
    #exit()
    C_points_table = np.load(path_1D + 'Cxq_points_table.npy', allow_pickle=True)
    print(C_points_table.shape)
    exit()
    point = np.array([alpha])
    C_lookup = interpn(C_points_table, C_values_table, point)
    return C_lookup

# print(Cx(alpha_in,beta_in,dh_in))
# print(Cy(alpha_in,beta_in))
#print(Cxq(alpha_in))







