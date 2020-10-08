#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 20:01:32 2020

@author: johnviljoen
"""

''' This script takes in the csv files extracted from the raw data in the form shown in the Nguyen
paper, and converts them into 3D, 2D, and 1D even grid data tables, filling the gaps filled in 
with linear interpolation '''

import numpy as np
import pandas as pd
from scipy.interpolate import griddata, interpn


path = '~/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen'

# In[import data]

############### aileron ################

aileron_path = path + '/raw_data/csv_processed/hifi_ailerons'

Cl_a20 = pd.read_csv(aileron_path + '/Cl_a20.csv', header=None)
Cl_a20_lef = pd.read_csv(aileron_path + '/Cl_a20_lef.csv', header=None)
Cn_a20 = pd.read_csv(aileron_path + '/Cn_a20.csv', header=None)
Cn_a20_lef = pd.read_csv(aileron_path + '/Cn_a20_lef.csv', header=None)
Cy_a20 = pd.read_csv(aileron_path + '/Cy_a20.csv', header=None)
Cy_a20_lef = pd.read_csv(aileron_path + '/Cy_a20_lef.csv', header=None)

############### C ################

C_path = path + '/raw_data/csv_processed/hifi_C'

Cl_dh_neg25 = pd.read_csv(C_path + '/Cl_dh_neg25.csv', header=None)
Cl_dh_0 = pd.read_csv(C_path + '/Cl_dh_0.csv', header=None)
Cl_dh_pos25 = pd.read_csv(C_path + '/Cl_dh_pos25.csv', header=None)

Cm_dh_neg25 = pd.read_csv(C_path + '/Cm_dh_neg25.csv', header=None)
Cm_dh_0 = pd.read_csv(C_path + '/Cm_dh_0.csv', header=None)
Cm_dh_pos25 = pd.read_csv(C_path + '/Cm_dh_pos25.csv', header=None)

Cn_dh_neg25 = pd.read_csv(C_path + '/Cn_dh_neg25.csv', header=None)
Cn_dh_0 = pd.read_csv(C_path + '/Cn_dh_0.csv', header=None)
Cn_dh_pos25 = pd.read_csv(C_path + '/Cn_dh_pos25.csv', header=None)

Cx_dh_neg25 = pd.read_csv(C_path + '/Cx_dh_neg25.csv', header=None)
Cx_dh_neg10 = pd.read_csv(C_path + '/Cx_dh_neg10.csv', header=None)
Cx_dh_0 = pd.read_csv(C_path + '/Cx_dh_0.csv', header=None)
Cx_dh_pos10 = pd.read_csv(C_path + '/Cx_dh_pos10.csv', header=None)
Cx_dh_pos25 = pd.read_csv(C_path + '/Cx_dh_pos25.csv', header=None)

Cy = pd.read_csv(C_path + '/Cy.csv', header=None)

Cz_dh_neg25 = pd.read_csv(C_path + '/Cz_dh_neg25.csv', header=None)
Cz_dh_neg10 = pd.read_csv(C_path + '/Cz_dh_neg10.csv', header=None)
Cz_dh_0 = pd.read_csv(C_path + '/Cz_dh_0.csv', header=None)
Cz_dh_pos10 = pd.read_csv(C_path + '/Cz_dh_pos10.csv', header=None)
Cz_dh_pos25 = pd.read_csv(C_path + '/Cz_dh_pos25.csv', header=None)

############### C_lef ################

C_lef_path = path + '/raw_data/csv_processed/hifi_C_lef'

Cx_lef = pd.read_csv(C_lef_path + '/Cx_lef.csv', header=None)
Cy_lef = pd.read_csv(C_lef_path + '/Cy_lef.csv', header=None)
Cz_lef = pd.read_csv(C_lef_path + '/Cz_lef.csv', header=None)
Cl_lef = pd.read_csv(C_lef_path + '/Cl_lef.csv', header=None)
Cm_lef = pd.read_csv(C_lef_path + '/Cm_lef.csv', header=None)
Cn_lef = pd.read_csv(C_lef_path + '/Cn_lef.csv', header=None)

############### damping ################

damping_path = path + '/raw_data/csv_processed/hifi_damping'

Clp = pd.read_csv(damping_path + '/Clp.csv', header=None)
Clr = pd.read_csv(damping_path + '/Clr.csv', header=None)
Cmq = pd.read_csv(damping_path + '/Cmq.csv', header=None)
Cnp = pd.read_csv(damping_path + '/Cnp.csv', header=None)
Cnr = pd.read_csv(damping_path + '/Cnr.csv', header=None)
Cxq = pd.read_csv(damping_path + '/Cxq.csv', header=None)
Cyp = pd.read_csv(damping_path + '/Cyp.csv', header=None)
Cyr = pd.read_csv(damping_path + '/Cyr.csv', header=None)
Czq = pd.read_csv(damping_path + '/Czq.csv', header=None)

############### damping_lef ################

damping_lef_path = path + '/raw_data/csv_processed/hifi_damping_lef'

delta_Clp_lef = pd.read_csv(damping_lef_path + '/delta_Clp_lef.csv', header=None)
delta_Clr_lef = pd.read_csv(damping_lef_path + '/delta_Clr_lef.csv', header=None)
delta_Cmq_lef = pd.read_csv(damping_lef_path + '/delta_Cmq_lef.csv', header=None)
delta_Cnp_lef = pd.read_csv(damping_lef_path + '/delta_Cnp_lef.csv', header=None)
delta_Cnr_lef = pd.read_csv(damping_lef_path + '/delta_Cnr_lef.csv', header=None)
delta_Cxq_lef = pd.read_csv(damping_lef_path + '/delta_Cxq_lef.csv', header=None)
delta_Cyp_lef = pd.read_csv(damping_lef_path + '/delta_Cyp_lef.csv', header=None)
delta_Cyr_lef = pd.read_csv(damping_lef_path + '/delta_Cyr_lef.csv', header=None)
delta_Czq_lef = pd.read_csv(damping_lef_path + '/delta_Czq_lef.csv', header=None)

############### other_coeffs ################

other_coeffs_path = path + '/raw_data/csv_processed/hifi_other_coeffs'

delta_Cl_beta = pd.read_csv(other_coeffs_path + '/delta_Cl_beta.csv', header=None)
delta_Cm = pd.read_csv(other_coeffs_path + '/delta_Cm.csv', header=None)
delta_Cn_beta = pd.read_csv(other_coeffs_path + '/delta_Cn_beta.csv', header=None)
eta_el = pd.read_csv(other_coeffs_path + '/eta_el.csv', header=None)

############### rudder ################

rudder_path = path + '/raw_data/csv_processed/hifi_rudder'

Cl_r30 = pd.read_csv(rudder_path + '/Cl_r30.csv', header=None)
Cn_r30 = pd.read_csv(rudder_path + '/Cn_r30.csv', header=None)
Cy_r30 = pd.read_csv(rudder_path + '/Cy_r30.csv', header=None)

# In[Unify 3D data]

########## for Cx ############
Cx_raw = np.zeros([20,19,5]) # initialise

Cx_raw[:,:,0] = Cx_dh_neg25
Cx_raw[:,:,1] = Cx_dh_neg10
Cx_raw[:,:,2] = Cx_dh_0
Cx_raw[:,:,3] = Cx_dh_pos10
Cx_raw[:,:,4] = Cx_dh_pos25

########## for Cz #############
Cz_raw = np.zeros([20,19,5]) # initialise

Cz_raw[:,:,0] = Cz_dh_neg25
Cz_raw[:,:,1] = Cz_dh_neg10
Cz_raw[:,:,2] = Cz_dh_0
Cz_raw[:,:,3] = Cz_dh_pos10
Cz_raw[:,:,4] = Cz_dh_pos25

########## for Cl #############
Cl_raw = np.zeros([20,19,3])

Cl_raw[:,:,0] = Cl_dh_neg25
Cl_raw[:,:,1] = Cl_dh_0
Cl_raw[:,:,2] = Cl_dh_pos25

########## for Cm #############
Cm_raw = np.zeros([20,19,3])

Cm_raw[:,:,0] = Cm_dh_neg25
Cm_raw[:,:,1] = Cm_dh_0
Cm_raw[:,:,2] = Cm_dh_pos25

########## for Cn #############
Cn_raw = np.zeros([20,19,3])

Cn_raw[:,:,0] = Cn_dh_neg25
Cn_raw[:,:,1] = Cn_dh_0
Cn_raw[:,:,2] = Cn_dh_pos25

# In[3D interpolate]

''' In this section the unequally spaced grid data is made into high fidelity
equally spaced data by using scipy.interpolate.griddata linear interpolation function
across the 3D datasets, this will allow for faster interpolation during simulation '''


# the C()_lef tables contain a smaller range of alpha than the other pieces of data
alphas_large_deg = np.array([-20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 90.0])
alphas_small_deg = np.array([-20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0])

# All tables utilise the same range of beta
beta_deg = np.array([-30.0, -25.0, -20.0, -15.0, -10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 25.0, 30.0])

# Some tables contain 5 values of dh, and some 3, and some 1
dh_large_deg = np.array([-25.0, -10.0, 0.0, 10.0, 25.0])
dh_small_deg = np.array([-25.0, 0.0, 25.0])

def process_3d_table_linear(alphas, betas, dhs, C):
    
    alphas_len = np.int(np.amax(alphas)-np.amin(alphas))+1
    alphas_min = np.amin(alphas)
    alphas_max = np.amax(alphas)
    alphas_even_deg = np.linspace(alphas_min,alphas_max,alphas_len)
    
    betas_len = np.int(np.amax(betas)-np.amin(betas))+1
    betas_min = np.amin(betas)
    betas_max = np.amax(betas)
    betas_even_deg = np.linspace(betas_min,betas_max,betas_len)
    
    dhs_len = np.int(np.amax(dhs)-np.amin(dhs))+1
    dhs_min = np.amin(dhs)
    dhs_max = np.amax(dhs)
    dhs_even_deg = np.linspace(dhs_min,dhs_max,dhs_len)
    
    C_processing = np.empty([alphas_len,betas_len,dhs_len]) # 3d scalar field (111 alphas, 61 betas, 5 dhs, 1 value)
    C_processing[:,:,:] = np.NaN
    
    points = np.empty([len(alphas),len(betas),len(dhs),3]) # initialise array for points we know values at for griddata
    points[:,:,:,:] = np.NaN
    
    values = np.empty([len(alphas),len(betas),len(dhs),1]) # initialise array for values of these points for griddata
    values[:,:,:,:] = np.NaN
    
    for i in range(0,alphas_len,1):
        print(i)
        #if alphas_even_deg[i] % 5. == 0: # checks to see if the current alphas value is a multiple of 5
        alpha_idx = np.where(alphas == alphas_even_deg[i])[0]
        # now that we have the idx at which the values occur at in the raw data
        # we can import that into the large Cx we are making
        #print(Cx_dh_neg25.values[idx])
        for j in range(0,betas_len,1):
            beta_idx = np.where(betas == betas_even_deg[j])[0]
            #print(beta_idx)
            #print(Cx_dh_neg25.values[:,beta_idx])
            for k in range(0,dhs_len,1):
                dh_idx = np.where(dhs == dhs_even_deg[k])[0]
                #print(dh_idx)
                
                if alpha_idx.size == 1 and beta_idx.size == 1 and dh_idx.size == 1: # we've found a point!
                    #dh_idx = k
                    C_processing[i,j,k] = C[alpha_idx,beta_idx,dh_idx]
                    points[alpha_idx,beta_idx,dh_idx,:] = [i-20,j-30,k-25]
                    values[alpha_idx,beta_idx,dh_idx,:] = C_processing[i,j,k]
    
    flattened_len = len(alphas)*len(betas)*len(dhs)
    points = np.reshape(points,[flattened_len,3]) # 3 dimensions
    values = np.reshape(values,[flattened_len,1]) # 1 dimension
    
    
    x = np.linspace(alphas_min,alphas_max,alphas_len)
    y = np.linspace(betas_min,betas_max,betas_len)
    z = np.linspace(dhs_min,dhs_max,dhs_len)
    xv, yv, zv = np.meshgrid(x,y,z)
    
    C_out = griddata(points, values, (xv,yv,zv), 'linear')
    #C_out = np.reshape(C_out,[alphas_len, betas_len, dhs_len,1]) # led to an incorrect reshuffling
    return C_out, (y, x, z)
    
Cx_values_table, Cx_points_table = process_3d_table_linear(alphas_large_deg, beta_deg, dh_large_deg, Cx_raw)
Cz_values_table, Cz_points_table = process_3d_table_linear(alphas_large_deg, beta_deg, dh_large_deg, Cz_raw)
Cl_values_table, Cl_points_table = process_3d_table_linear(alphas_large_deg, beta_deg, dh_small_deg, Cl_raw)
Cm_values_table, Cm_points_table = process_3d_table_linear(alphas_large_deg, beta_deg, dh_small_deg, Cm_raw)
Cn_values_table, Cn_points_table = process_3d_table_linear(alphas_large_deg, beta_deg, dh_small_deg, Cn_raw)

# In[2D interpolate]

''' In this section the unequally spaced grid data is made into high fidelity
equally spaced data by using scipy.interpolate.griddata linear interpolation function
across the 2D datasets, this will allow for faster interpolation during simulation '''

def process_2d_table_linear(alphas, betas, C):
    
    alphas_len = np.int(np.amax(alphas)-np.amin(alphas))+1
    alphas_min = np.amin(alphas)
    alphas_max = np.amax(alphas)
    alphas_even_deg = np.linspace(alphas_min,alphas_max,alphas_len)
    
    betas_len = np.int(np.amax(betas)-np.amin(betas))+1
    betas_min = np.amin(betas)
    betas_max = np.amax(betas)
    betas_even_deg = np.linspace(betas_min,betas_max,betas_len)
    
    C_processing = np.empty([alphas_len,betas_len]) # 3d scalar field (111 alphas, 61 betas, 5 dhs, 1 value)
    C_processing[:,:] = np.NaN
    
    points = np.empty([len(alphas),len(betas),2]) # initialise array for points we know values at for griddata
    points[:,:,:] = np.NaN
    
    values = np.empty([len(alphas),len(betas),1]) # initialise array for values of these points for griddata
    values[:,:,:] = np.NaN
    
    for i in range(0,alphas_len,1):
        print(i)
        #if alphas_even_deg[i] % 5. == 0: # checks to see if the current alphas value is a multiple of 5
        alpha_idx = np.where(alphas == alphas_even_deg[i])[0]
        # now that we have the idx at which the values occur at in the raw data
        # we can import that into the large Cx we are making
        #print(Cx_dh_neg25.values[idx])
        for j in range(0,betas_len,1):
            beta_idx = np.where(betas == betas_even_deg[j])[0]
            #print(beta_idx)
            #print(Cx_dh_neg25.values[:,beta_idx])
            if alpha_idx.size == 1 and beta_idx.size == 1: # we've found a point!
                #dh_idx = k
                C_processing[i,j] = C.to_numpy()[alpha_idx,beta_idx]
                points[alpha_idx,beta_idx,:] = [i-20,j-30]
                values[alpha_idx,beta_idx,:] = C_processing[i,j]
    
    flattened_len = len(alphas)*len(betas)
    points = np.reshape(points,[flattened_len,2]) # 2 dimensions
    values = np.reshape(values,[flattened_len,1]) # 1 dimension
    
    
    x = np.linspace(alphas_min,alphas_max,alphas_len)
    y = np.linspace(betas_min,betas_max,betas_len)
    xv, yv = np.meshgrid(x,y)
    
    C_out = griddata(points, values, (xv,yv), 'linear')
    #C_out = np.reshape(C_out,[alphas_len, betas_len, dhs_len,1]) # led to an incorrect reshuffling
    return C_out, (y, x)    

#Cy
Cy_values_table, Cy_points_table = process_2d_table_linear(alphas_large_deg,beta_deg,Cy)

#C()lef's
Cx_lef_values_table, Cx_lef_points_table = process_2d_table_linear(alphas_small_deg,beta_deg,Cx_lef)
Cy_lef_values_table, Cy_lef_points_table = process_2d_table_linear(alphas_small_deg,beta_deg,Cy_lef)
Cz_lef_values_table, Cz_lef_points_table = process_2d_table_linear(alphas_small_deg,beta_deg,Cz_lef)
Cl_lef_values_table, Cl_lef_points_table = process_2d_table_linear(alphas_small_deg,beta_deg,Cl_lef)
Cm_lef_values_table, Cm_lef_points_table = process_2d_table_linear(alphas_small_deg,beta_deg,Cm_lef)
Cn_lef_values_table, Cn_lef_points_table = process_2d_table_linear(alphas_small_deg,beta_deg,Cn_lef)

#C()_a20's
Cy_a20_values_table, Cy_a20_points_table = process_2d_table_linear(alphas_large_deg,beta_deg,Cy_a20)
Cl_a20_values_table, Cl_a20_points_table = process_2d_table_linear(alphas_large_deg,beta_deg,Cl_a20)
Cn_a20_values_table, Cn_a20_points_table = process_2d_table_linear(alphas_large_deg,beta_deg,Cn_a20)

#C()_a20_lef's
Cy_a20_lef_values_table, Cy_a20_lef_points_table = process_2d_table_linear(alphas_small_deg,beta_deg,Cy_a20_lef)
Cl_a20_lef_values_table, Cl_a20_lef_points_table = process_2d_table_linear(alphas_small_deg,beta_deg,Cl_a20_lef)
Cn_a20_lef_values_table, Cn_a20_lef_points_table = process_2d_table_linear(alphas_small_deg,beta_deg,Cn_a20_lef)

#C()_r30's
Cy_r30_values_table, Cy_r30_points_table = process_2d_table_linear(alphas_large_deg,beta_deg,Cy_r30)
Cl_r30_values_table, Cl_r30_points_table = process_2d_table_linear(alphas_large_deg,beta_deg,Cl_r30)
Cn_r30_values_table, Cn_r30_points_table = process_2d_table_linear(alphas_large_deg,beta_deg,Cn_r30)

# In[1D interpolate]
    
''' In this section the unequally spaced grid data is made into high fidelity
equally spaced data by using scipy.interpolate.griddata linear interpolation function
across the 1D datasets, this will allow for faster interpolation during simulation '''

def process_1d_table_linear(alphas, C):
    
    alphas_len = np.int(np.amax(alphas)-np.amin(alphas))+1
    alphas_min = np.amin(alphas)
    alphas_max = np.amax(alphas)
    alphas_even_deg = np.linspace(alphas_min,alphas_max,alphas_len)
    
    C_processing = np.empty([alphas_len]) # 3d scalar field (111 alphas, 61 betas, 5 dhs, 1 value)
    C_processing[:] = np.NaN
    
    points = np.empty([len(alphas),1]) # initialise array for points we know values at for griddata
    points[:,:] = np.NaN
    
    values = np.empty([len(alphas),1]) # initialise array for values of these points for griddata
    values[:,:] = np.NaN
    
    for i in range(0,alphas_len,1):
        print(i)
        #if alphas_even_deg[i] % 5. == 0: # checks to see if the current alphas value is a multiple of 5
        alpha_idx = np.where(alphas == alphas_even_deg[i])[0]
        # now that we have the idx at which the values occur at in the raw data
        # we can import that into the large Cx we are making
        #print(Cx_dh_neg25.values[idx])
        if alpha_idx.size == 1: # we've found a point!
            #dh_idx = k
            C_processing[i] = C.to_numpy()[alpha_idx]
            points[alpha_idx,:] = [i-20]
            values[alpha_idx,:] = C_processing[i]
    
    flattened_len = len(alphas)
    points = np.reshape(points,[flattened_len,1]) # 1 dimensions
    values = np.reshape(values,[flattened_len,1]) # 1 dimension
    
    
    x = np.linspace(alphas_min,alphas_max,alphas_len)
    xv = np.meshgrid(x)
    
    C_out = griddata(points, values, (xv), 'linear')
    #C_out = np.reshape(C_out,[alphas_len, betas_len, dhs_len,1]) # led to an incorrect reshuffling
    return C_out, (x)    
    
# p derivatives
Cyp_values_table, Cyp_points_table = process_1d_table_linear(alphas_large_deg,Cyp)
Cnp_values_table, Cnp_points_table = process_1d_table_linear(alphas_large_deg,Cnp)
Clp_values_table, Clp_points_table = process_1d_table_linear(alphas_large_deg,Clp)

delta_Cyp_lef_values_table, delta_Cyp_lef_points_table = process_1d_table_linear(alphas_small_deg,delta_Cyp_lef)
delta_Cnp_lef_values_table, delta_Cnp_lef_points_table = process_1d_table_linear(alphas_small_deg,delta_Cnp_lef)
delta_Clp_lef_values_table, delta_Clp_lef_points_table = process_1d_table_linear(alphas_small_deg,delta_Clp_lef)

# q derivatives
Cxq_values_table, Cxq_points_table = process_1d_table_linear(alphas_large_deg,Cxq)
Czq_values_table, Czq_points_table = process_1d_table_linear(alphas_large_deg,Czq)
Cmq_values_table, Cmq_points_table = process_1d_table_linear(alphas_large_deg,Cmq)

delta_Cxq_lef_values_table, delta_Cxq_lef_points_table = process_1d_table_linear(alphas_small_deg,delta_Cxq_lef)
delta_Czq_lef_values_table, delta_Czq_lef_points_table = process_1d_table_linear(alphas_small_deg,delta_Czq_lef)
delta_Cmq_lef_values_table, delta_Cmq_lef_points_table = process_1d_table_linear(alphas_small_deg,delta_Cmq_lef)

# r derivatives
Cyr_values_table, Cyr_points_table = process_1d_table_linear(alphas_large_deg,Cyr)
Cnr_values_table, Cnr_points_table = process_1d_table_linear(alphas_large_deg,Cnr)
Clr_values_table, Clr_points_table = process_1d_table_linear(alphas_large_deg,Clr)

delta_Cyr_lef_values_table, delta_Cyr_lef_points_table = process_1d_table_linear(alphas_small_deg,delta_Cyr_lef)
delta_Cnr_lef_values_table, delta_Cnr_lef_points_table = process_1d_table_linear(alphas_small_deg,delta_Cnr_lef)
delta_Clr_lef_values_table, delta_Clr_lef_points_table = process_1d_table_linear(alphas_small_deg,delta_Clr_lef)

# other
delta_Cm_values_table, delta_Cm_points_table = process_1d_table_linear(alphas_large_deg,delta_Cm)
delta_Cn_beta_values_table, delta_Cn_beta_points_table = process_1d_table_linear(alphas_large_deg,delta_Cn_beta)
#delta_Cn_da_values_table, delta_Cn_da_points_table = process_1d_table_linear(alphas_large_deg,delta_Cn_da)
delta_Cl_beta_values_table, delta_Cl_beta_points_table = process_1d_table_linear(alphas_large_deg,delta_Cl_beta)

# In[save processed data]

################# 3D ##################

path_3D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/3D/'

np.save(path_3D + 'Cx_values_table.npy',Cx_values_table)
np.save(path_3D + 'Cz_values_table.npy',Cz_values_table)
np.save(path_3D + 'Cl_values_table.npy',Cl_values_table)
np.save(path_3D + 'Cm_values_table.npy',Cm_values_table)
np.save(path_3D + 'Cn_values_table.npy',Cn_values_table)
np.save(path_3D + 'Cx_points_table.npy',Cx_points_table)
np.save(path_3D + 'Cz_points_table.npy',Cz_points_table)
np.save(path_3D + 'Cl_points_table.npy',Cl_points_table)
np.save(path_3D + 'Cm_points_table.npy',Cm_points_table)
np.save(path_3D + 'Cn_points_table.npy',Cn_points_table)

################## 2D ##################

path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'

# Cy
np.save(path_2D + 'Cy_values_table.npy',Cy_values_table)
np.save(path_2D + 'Cy_points_table.npy',Cy_points_table)

#C()lef's
np.save(path_2D + 'Cx_lef_values_table.npy',Cx_lef_values_table)
np.save(path_2D + 'Cx_lef_points_table.npy',Cx_lef_points_table)
np.save(path_2D + 'Cy_lef_values_table.npy',Cy_lef_values_table)
np.save(path_2D + 'Cy_lef_points_table.npy',Cy_lef_points_table)
np.save(path_2D + 'Cz_lef_values_table.npy',Cz_lef_values_table)
np.save(path_2D + 'Cz_lef_points_table.npy',Cz_lef_points_table)
np.save(path_2D + 'Cl_lef_values_table.npy',Cl_lef_values_table)
np.save(path_2D + 'Cl_lef_points_table.npy',Cl_lef_points_table)
np.save(path_2D + 'Cm_lef_values_table.npy',Cm_lef_values_table)
np.save(path_2D + 'Cm_lef_points_table.npy',Cm_lef_points_table)
np.save(path_2D + 'Cn_lef_values_table.npy',Cn_lef_values_table)
np.save(path_2D + 'Cn_lef_points_table.npy',Cn_lef_points_table)

#C()_a20's
np.save(path_2D + 'Cy_a20_values_table.npy',Cy_a20_values_table)
np.save(path_2D + 'Cy_a20_points_table.npy',Cy_a20_points_table)
np.save(path_2D + 'Cl_a20_values_table.npy',Cl_a20_values_table)
np.save(path_2D + 'Cl_a20_points_table.npy',Cl_a20_points_table)
np.save(path_2D + 'Cn_a20_values_table.npy',Cn_a20_values_table)
np.save(path_2D + 'Cn_a20_points_table.npy',Cn_a20_points_table)

#C()_a20_lef's
np.save(path_2D + 'Cy_a20_lef_values_table.npy',Cy_a20_lef_values_table)
np.save(path_2D + 'Cy_a20_lef_points_table.npy',Cy_a20_lef_points_table)
np.save(path_2D + 'Cl_a20_lef_values_table.npy',Cl_a20_lef_values_table)
np.save(path_2D + 'Cl_a20_lef_points_table.npy',Cl_a20_lef_points_table)
np.save(path_2D + 'Cn_a20_lef_values_table.npy',Cn_a20_lef_values_table)
np.save(path_2D + 'Cn_a20_lef_points_table.npy',Cn_a20_lef_points_table)

#C()_r30's
np.save(path_2D + 'Cy_r30_values_table.npy',Cy_r30_values_table)
np.save(path_2D + 'Cy_r30_points_table.npy',Cy_r30_points_table)
np.save(path_2D + 'Cl_r30_values_table.npy',Cl_r30_values_table)
np.save(path_2D + 'Cl_r30_points_table.npy',Cl_r30_points_table)
np.save(path_2D + 'Cn_r30_values_table.npy',Cn_r30_values_table)
np.save(path_2D + 'Cn_r30_points_table.npy',Cn_r30_points_table)

#################### 1D ####################

path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'

# p derivatives
np.save(path_1D + 'Cyp_values_table.npy',Cyp_values_table)
np.save(path_1D + 'Cyp_points_table.npy',Cyp_points_table)
np.save(path_1D + 'Cnp_values_table.npy',Cnp_values_table)
np.save(path_1D + 'Cnp_points_table.npy',Cnp_points_table)
np.save(path_1D + 'Clp_values_table.npy',Clp_values_table)
np.save(path_1D + 'Clp_points_table.npy',Clp_points_table)

np.save(path_1D + 'delta_Cyp_lef_values_table.npy',delta_Cyp_lef_values_table)
np.save(path_1D + 'delta_Cyp_lef_points_table.npy',delta_Cyp_lef_points_table)
np.save(path_1D + 'delta_Cyp_lef_values_table.npy',delta_Cyp_lef_values_table)
np.save(path_1D + 'delta_Cyp_lef_points_table.npy',delta_Cyp_lef_points_table)
np.save(path_1D + 'delta_Cyp_lef_values_table.npy',delta_Cyp_lef_values_table)
np.save(path_1D + 'delta_Cyp_lef_points_table.npy',delta_Cyp_lef_points_table)

# q derivatives
np.save(path_1D + 'Cxq_values_table.npy',Cxq_values_table)
np.save(path_1D + 'Cxq_points_table.npy',Cxq_points_table)
np.save(path_1D + 'Czq_values_table.npy',Czq_values_table)
np.save(path_1D + 'Czq_points_table.npy',Czq_points_table)
np.save(path_1D + 'Cmq_values_table.npy',Cmq_values_table)
np.save(path_1D + 'Cmq_points_table.npy',Cmq_points_table)

np.save(path_1D + 'delta_Cxq_lef_values_table.npy',delta_Cxq_lef_values_table)
np.save(path_1D + 'delta_Cxq_lef_points_table.npy',delta_Cxq_lef_points_table)
np.save(path_1D + 'delta_Czq_lef_values_table.npy',delta_Czq_lef_values_table)
np.save(path_1D + 'delta_Czq_lef_points_table.npy',delta_Czq_lef_points_table)
np.save(path_1D + 'delta_Cmq_lef_values_table.npy',delta_Cmq_lef_values_table)
np.save(path_1D + 'delta_Cmq_lef_points_table.npy',delta_Cmq_lef_points_table)

# r derivatives
np.save(path_1D + 'Cyr_values_table.npy',Cyr_values_table)
np.save(path_1D + 'Cyr_points_table.npy',Cyr_points_table)
np.save(path_1D + 'Cnr_values_table.npy',Cnr_values_table)
np.save(path_1D + 'Cnr_points_table.npy',Cnr_points_table)
np.save(path_1D + 'Clr_values_table.npy',Clr_values_table)
np.save(path_1D + 'Clr_points_table.npy',Clr_points_table)

np.save(path_1D + 'delta_Cyr_lef_values_table.npy',delta_Cyr_lef_values_table)
np.save(path_1D + 'delta_Cyr_lef_points_table.npy',delta_Cyr_lef_points_table)
np.save(path_1D + 'delta_Cnr_lef_values_table.npy',delta_Cnr_lef_values_table)
np.save(path_1D + 'delta_Cnr_lef_points_table.npy',delta_Cnr_lef_points_table)
np.save(path_1D + 'delta_Clr_lef_values_table.npy',delta_Clr_lef_values_table)
np.save(path_1D + 'delta_Clr_lef_points_table.npy',delta_Clr_lef_points_table)

# other
np.save(path_1D + 'delta_Cm_values_table.npy',delta_Cm_values_table)
np.save(path_1D + 'delta_Cm_points_table.npy',delta_Cm_points_table)
np.save(path_1D + 'delta_Cn_beta_values_table.npy',delta_Cn_beta_values_table)
np.save(path_1D + 'delta_Cn_beta_points_table.npy',delta_Cn_beta_points_table)
np.save(path_1D + 'delta_Cl_beta_values_table.npy',delta_Cl_beta_values_table)
np.save(path_1D + 'delta_Cl_beta_points_table.npy',delta_Cl_beta_points_table)

# In[]test interpn on this grid
alpha_in = -19.23
beta_in = -28.42
dh_in = -23.2

point = np.array([beta_in,alpha_in,dh_in])

print(interpn(Cx_points_table, Cx_values_table[:,:,:,0], point))
print(interpn(Cz_points_table, Cz_values_table[:,:,:,0], point))
print(interpn(Cl_points_table, Cl_values_table[:,:,:,0], point))
print(interpn(Cm_points_table, Cm_values_table[:,:,:,0], point))
print(interpn(Cn_points_table, Cn_values_table[:,:,:,0], point))

def Cx(alpha,beta,dh):
    point = np.array([beta, alpha, dh])
    Cx_lookup = interpn(Cx_points_table, Cx_values_table[:,:,:,0], point)
    return Cx_lookup

def Cz(alpha,beta,dh):
    point = np.array([beta, alpha, dh])
    Cz_lookup = interpn(Cz_points_table, Cz_values_table[:,:,:,0], point)
    return Cz_lookup

print(Cx(alpha_in,beta_in,dh_in))
print(Cz(alpha_in,beta_in,dh_in))


