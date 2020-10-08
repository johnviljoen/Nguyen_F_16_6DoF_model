#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 17:39:34 2020

@author: johnviljoen
"""

# In[dictionary of loading aerodynamic data]

import numpy as np

################# 3D ##################

path_3D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/3D/'

Cx_values_table = np.load(path_3D + 'Cx_values_table.npy')
Cz_values_table = np.load(path_3D + 'Cz_values_table.npy')
Cl_values_table = np.load(path_3D + 'Cl_values_table.npy')
Cm_values_table = np.load(path_3D + 'Cm_values_table.npy')
Cn_values_table = np.load(path_3D + 'Cn_values_table.npy')
Cx_points_table = np.load(path_3D + 'Cx_points_table.npy', allow_pickle=True)
Cz_points_table = np.load(path_3D + 'Cz_points_table.npy', allow_pickle=True)
Cl_points_table = np.load(path_3D + 'Cl_points_table.npy', allow_pickle=True)
Cm_points_table = np.load(path_3D + 'Cm_points_table.npy', allow_pickle=True)
Cn_points_table = np.load(path_3D + 'Cn_points_table.npy', allow_pickle=True)

################## 2D ##################

path_2D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/2D/'

# Cy
Cy_values_table = np.load(path_2D + 'Cy_values_table.npy')
Cy_points_table = np.load(path_2D + 'Cy_points_table.npy', allow_pickle=True)

#C()lef's
Cx_lef_values_table = np.load(path_2D + 'Cx_lef_values_table.npy')
Cx_lef_points_table = np.load(path_2D + 'Cx_lef_points_table.npy', allow_pickle=True)
Cy_lef_values_table = np.load(path_2D + 'Cy_lef_values_table.npy')
Cy_lef_points_table = np.load(path_2D + 'Cy_lef_points_table.npy', allow_pickle=True)
Cz_lef_values_table = np.load(path_2D + 'Cz_lef_values_table.npy')
Cz_lef_points_table = np.load(path_2D + 'Cz_lef_points_table.npy', allow_pickle=True)
Cl_lef_values_table = np.load(path_2D + 'Cl_lef_values_table.npy')
Cl_lef_points_table = np.load(path_2D + 'Cl_lef_points_table.npy', allow_pickle=True)
Cm_lef_values_table = np.load(path_2D + 'Cm_lef_values_table.npy')
Cm_lef_points_table = np.load(path_2D + 'Cm_lef_points_table.npy', allow_pickle=True)
Cn_lef_values_table = np.load(path_2D + 'Cn_lef_values_table.npy')
Cn_lef_points_table = np.load(path_2D + 'Cn_lef_points_table.npy', allow_pickle=True)

#C()_a20's
Cy_a20_values_table = np.load(path_2D + 'Cy_a20_values_table.npy')
Cy_a20_points_table = np.load(path_2D + 'Cy_a20_points_table.npy', allow_pickle=True)
Cl_a20_values_table = np.load(path_2D + 'Cl_a20_values_table.npy')
Cl_a20_points_table = np.load(path_2D + 'Cl_a20_points_table.npy', allow_pickle=True)
Cn_a20_values_table = np.load(path_2D + 'Cn_a20_values_table.npy')
Cn_a20_points_table = np.load(path_2D + 'Cn_a20_points_table.npy', allow_pickle=True)

#C()_a20_lef's
Cy_a20_lef_values_table = np.load(path_2D + 'Cy_a20_lef_values_table.npy')
Cy_a20_lef_points_table = np.load(path_2D + 'Cy_a20_lef_points_table.npy', allow_pickle=True)
Cl_a20_lef_values_table = np.load(path_2D + 'Cl_a20_lef_values_table.npy')
Cl_a20_lef_points_table = np.load(path_2D + 'Cl_a20_lef_points_table.npy', allow_pickle=True)
Cn_a20_lef_values_table = np.load(path_2D + 'Cn_a20_lef_values_table.npy')
Cn_a20_lef_points_table = np.load(path_2D + 'Cn_a20_lef_points_table.npy', allow_pickle=True)

#C()_r30's
Cy_r30_values_table = np.load(path_2D + 'Cy_r30_values_table.npy')
Cy_r30_points_table = np.load(path_2D + 'Cy_r30_points_table.npy', allow_pickle=True)
Cl_r30_values_table = np.load(path_2D + 'Cl_r30_values_table.npy')
Cl_r30_points_table = np.load(path_2D + 'Cl_r30_points_table.npy', allow_pickle=True)
Cn_r30_values_table = np.load(path_2D + 'Cn_r30_values_table.npy')
Cn_r30_points_table = np.load(path_2D + 'Cn_r30_points_table.npy', allow_pickle=True)

#################### 1D ####################

path_1D = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/processed_data/1D/'

# p derivatives
Cyp_values_table = np.load(path_1D + 'Cyp_values_table.npy')
Cyp_points_table = np.load(path_1D + 'Cyp_points_table.npy', allow_pickle=True)
Cnp_values_table = np.load(path_1D + 'Cnp_values_table.npy')
Cnp_points_table = np.load(path_1D + 'Cnp_points_table.npy', allow_pickle=True)
Clp_values_table = np.load(path_1D + 'Clp_values_table.npy')
Clp_points_table = np.load(path_1D + 'Clp_points_table.npy', allow_pickle=True)

delta_Cyp_lef_values_table = np.load(path_1D + 'delta_Cyp_lef_values_table.npy')
delta_Cyp_lef_points_table = np.load(path_1D + 'delta_Cyp_lef_points_table.npy', allow_pickle=True)
delta_Cyp_lef_values_table = np.load(path_1D + 'delta_Cyp_lef_values_table.npy')
delta_Cyp_lef_points_table = np.load(path_1D + 'delta_Cyp_lef_points_table.npy', allow_pickle=True)
delta_Cyp_lef_values_table = np.load(path_1D + 'delta_Cyp_lef_values_table.npy')
delta_Cyp_lef_points_table = np.load(path_1D + 'delta_Cyp_lef_points_table.npy', allow_pickle=True)

# q derivatives
Cxq_values_table = np.load(path_1D + 'Cxq_values_table.npy')
Cxq_points_table = np.load(path_1D + 'Cxq_points_table.npy', allow_pickle=True)
Czq_values_table = np.load(path_1D + 'Czq_values_table.npy')
Czq_points_table = np.load(path_1D + 'Czq_points_table.npy', allow_pickle=True)
Cmq_values_table = np.load(path_1D + 'Cmq_values_table.npy')
Cmq_points_table = np.load(path_1D + 'Cmq_points_table.npy', allow_pickle=True)

delta_Cxq_lef_values_table = np.load(path_1D + 'delta_Cxq_lef_values_table.npy')
delta_Cxq_lef_points_table = np.load(path_1D + 'delta_Cxq_lef_points_table.npy', allow_pickle=True)
delta_Czq_lef_values_table = np.load(path_1D + 'delta_Czq_lef_values_table.npy')
delta_Czq_lef_points_table = np.load(path_1D + 'delta_Czq_lef_points_table.npy', allow_pickle=True)
delta_Cmq_lef_values_table = np.load(path_1D + 'delta_Cmq_lef_values_table.npy')
delta_Cmq_lef_points_table = np.load(path_1D + 'delta_Cmq_lef_points_table.npy', allow_pickle=True)

# r derivatives
Cyr_values_table = np.load(path_1D + 'Cyr_values_table.npy')
Cyr_points_table = np.load(path_1D + 'Cyr_points_table.npy', allow_pickle=True)
Cnr_values_table = np.load(path_1D + 'Cnr_values_table.npy')
Cnr_points_table = np.load(path_1D + 'Cnr_points_table.npy', allow_pickle=True)
Clr_values_table = np.load(path_1D + 'Clr_values_table.npy')
Clr_points_table = np.load(path_1D + 'Clr_points_table.npy', allow_pickle=True)

delta_Cyr_lef_values_table = np.load(path_1D + 'delta_Cyr_lef_values_table.npy')
delta_Cyr_lef_points_table = np.load(path_1D + 'delta_Cyr_lef_points_table.npy', allow_pickle=True)
delta_Cnr_lef_values_table = np.load(path_1D + 'delta_Cnr_lef_values_table.npy')
delta_Cnr_lef_points_table = np.load(path_1D + 'delta_Cnr_lef_points_table.npy', allow_pickle=True)
delta_Clr_lef_values_table = np.load(path_1D + 'delta_Clr_lef_values_table.npy')
delta_Clr_lef_points_table = np.load(path_1D + 'delta_Clr_lef_points_table.npy', allow_pickle=True)

# other
delta_Cm_values_table = np.load(path_1D + 'delta_Cm_values_table.npy')
delta_Cm_points_table = np.load(path_1D + 'delta_Cm_points_table.npy', allow_pickle=True)
delta_Cn_beta_values_table = np.load(path_1D + 'delta_Cn_beta_values_table.npy')
delta_Cn_beta_points_table = np.load(path_1D + 'delta_Cn_beta_points_table.npy', allow_pickle=True)
delta_Cl_beta_values_table = np.load(path_1D + 'delta_Cl_beta_values_table.npy')
delta_Cl_beta_points_table = np.load(path_1D + 'delta_Cl_beta_points_table.npy', allow_pickle=True)

# In[testing]

from scipy.interpolate import interp1d

points = Cyp_points_table
values = Cyp_values_table

values = np.reshape(values,[111])

f = interp1d(points, values)

point = (5)

print(f(point))