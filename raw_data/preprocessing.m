clear all

%%%%%%%%%%%%% ADD /raw TO THE PATH %%%%%%%%%%%%%%

write = 1; %
csv_path = '/home/johnviljoen/Documents/Aero Year 4/RP4/Python/Modelling/mk4 Nguyen/raw_data/csv_processed/';

%% CL
load('CL0120_ALPHA1_BETA1_DH2_601.txt')
Cl_old = CL0120_ALPHA1_BETA1_DH2_601;
Cl_dh_neg25 = reshape(Cl_old(1:380),[20,19]);
Cl_dh_0 = reshape(Cl_old(381:760),[20,19]);
Cl_dh_pos25 = reshape(Cl_old(761:1140),[20,19]);

%% Cm
load('CM0120_ALPHA1_BETA1_DH1_101.dat')
Cm_old = CM0120_ALPHA1_BETA1_DH1_101;
Cm_dh_neg25 = reshape(Cm_old(1:380),[20,19]);
Cm_dh_neg10 = reshape(Cm_old(381:760),[20,19]);
Cm_dh_0 = reshape(Cm_old(761:1140),[20,19]);
Cm_dh_pos10 = reshape(Cm_old(1141:1520),[20,19]);
Cm_dh_pos25 = reshape(Cm_old(1521:1900),[20,19]);

%% Cn
load('CN0120_ALPHA1_BETA1_DH2_501.dat')
Cn_old = CN0120_ALPHA1_BETA1_DH2_501;
Cn_dh_neg25 = reshape(Cn_old(1:380),[20,19]);
Cn_dh_0 = reshape(Cn_old(381:760),[20,19]);
Cn_dh_pos25 = reshape(Cn_old(761:1140),[20,19]);

%% Cx
load('CX0120_ALPHA1_BETA1_DH1_201.dat')
Cx_old = CX0120_ALPHA1_BETA1_DH1_201;
Cx_dh_neg25 = reshape(Cx_old(1:380),[20,19]);
Cx_dh_neg10 = reshape(Cx_old(381:760),[20,19]);
Cx_dh_0 = reshape(Cx_old(761:1140),[20,19]);
Cx_dh_pos10 = reshape(Cx_old(1141:1520),[20,19]);
Cx_dh_pos25 = reshape(Cx_old(1521:1900),[20,19]);

%% Cy
load('CY0320_ALPHA1_BETA1_401.dat')
Cy_old = CY0320_ALPHA1_BETA1_401; % the same one as was chosen by the C code
Cy = reshape(Cy_old,[20,19]);

%% Cz
load('CZ0120_ALPHA1_BETA1_DH1_301.dat')
Cz_old = CZ0120_ALPHA1_BETA1_DH1_301;
Cz_dh_neg25 = reshape(Cz_old(1:380),[20,19]);
Cz_dh_neg10 = reshape(Cz_old(381:760),[20,19]);
Cz_dh_0 = reshape(Cz_old(761:1140),[20,19]);
Cz_dh_pos10 = reshape(Cz_old(1141:1520),[20,19]);
Cz_dh_pos25 = reshape(Cz_old(1521:1900),[20,19]);

%% Cx_lef
load('CX0820_ALPHA2_BETA1_202.dat')
Cx_lef_old = CX0820_ALPHA2_BETA1_202;
Cx_lef = reshape(Cx_lef_old,[14,19]);

%% Cz_lef
load('CZ0820_ALPHA2_BETA1_302.dat')
Cz_lef_old = CZ0820_ALPHA2_BETA1_302;
Cz_lef = reshape(Cz_lef_old,[14,19]);

%% Cm_lef
load('CM0820_ALPHA2_BETA1_102.dat')
Cm_lef_old = CM0820_ALPHA2_BETA1_102;
Cm_lef = reshape(Cm_lef_old,[14,19]);

%% Cy_lef
load('CY0820_ALPHA2_BETA1_402.dat')
Cy_lef_old = CY0820_ALPHA2_BETA1_402;
Cy_lef = reshape(Cy_lef_old,[14,19]);

%% Cn_lef
load('CN0820_ALPHA2_BETA1_502.dat')
Cn_lef_old = CN0820_ALPHA2_BETA1_502;
Cn_lef = reshape(Cn_lef_old,[14,19]);

%% Cl_lef
load('CL0820_ALPHA2_BETA1_602.dat')
Cl_lef_old = CL0820_ALPHA2_BETA1_602;
Cl_lef = reshape(Cl_lef_old,[14,19]);

%% Cxq
load('CX1120_ALPHA1_204.dat')
Cxq_old = CX1120_ALPHA1_204;
Cxq = reshape(Cxq_old,[20,1]);

%% Czq
load('CZ1120_ALPHA1_304.dat')
Czq_old = CZ1120_ALPHA1_304;
Czq = reshape(Czq_old,[20,1]);

%% Cmq
load('CM1120_ALPHA1_104.dat')
Cmq_old = CM1120_ALPHA1_104;
Cmq = reshape(Cmq_old,[20,1]);

%% Cyp
load('CY1220_ALPHA1_408.dat')
Cyp_old = CY1220_ALPHA1_408;
Cyp = reshape(Cyp_old,[20,1]);

%% Cyr
load('CY1320_ALPHA1_406.dat')
Cyr_old = CY1320_ALPHA1_406;
Cyr = reshape(Cyr_old,[20,1]);

%% Cnr
load('CN1320_ALPHA1_506.dat')
Cnr_old = CN1320_ALPHA1_506;
Cnr = reshape(Cnr_old,[20,1]);

%% Cnp
load('CN1220_ALPHA1_508.dat')
Cnp_old = CN1220_ALPHA1_508;
Cnp = reshape(Cnp_old,[20,1]);

%% Clp
load('CL1220_ALPHA1_608.dat')
Clp_old = CL1220_ALPHA1_608;
Clp = reshape(Clp_old,[20,1]);

%% Clr
load('CL1320_ALPHA1_606.dat')
Clr_old = CL1320_ALPHA1_606;
Clr = reshape(Clr_old,[20,1]);

%% delta_Cxq_lef
load('CX1420_ALPHA2_205.dat')
delta_Cxq_lef_old = CX1420_ALPHA2_205;
delta_Cxq_lef = reshape(delta_Cxq_lef_old,[14,1]);

%% delta_Cyr_lef
load('CY1620_ALPHA2_407.dat')
delta_Cyr_lef_old = CY1620_ALPHA2_407;
delta_Cyr_lef = reshape(delta_Cyr_lef_old,[14,1]);

%% delta_Cyp_lef
load('CY1520_ALPHA2_409.dat')
delta_Cyp_lef_old = CY1520_ALPHA2_409;
delta_Cyp_lef = reshape(delta_Cyp_lef_old,[14,1]);

%% delta_Czq_lef
load('CZ1420_ALPHA2_305.dat')
delta_Czq_lef_old = CZ1420_ALPHA2_305;
delta_Czq_lef = reshape(delta_Czq_lef_old,[14,1]);

%% delta_Clr_lef
load('CL1620_ALPHA2_607.dat')
delta_Clr_lef_old = CL1620_ALPHA2_607;
delta_Clr_lef = reshape(delta_Clr_lef_old,[14,1]);

%% delta_Clp_lef
load('CL1520_ALPHA2_609.dat')
delta_Clp_lef_old = CL1520_ALPHA2_609;
delta_Clp_lef = reshape(delta_Clp_lef_old,[14,1]);

%% delta_Cmq_lef
load('CM1420_ALPHA2_105.dat')
delta_Cmq_lef_old = CM1420_ALPHA2_105;
delta_Cmq_lef = reshape(delta_Cmq_lef_old,[14,1]);

%% delta_Cnr_lef
load('CN1620_ALPHA2_507.dat')
delta_Cnr_lef_old = CN1620_ALPHA2_507;
delta_Cnr_lef = reshape(delta_Cnr_lef_old,[14,1]);

%% delta_Cnp_lef
load('CN1520_ALPHA2_509.dat')
delta_Cnp_lef_old = CN1520_ALPHA2_509;
delta_Cnp_lef = reshape(delta_Cnp_lef_old,[14,1]);

%% Cy_r30
load('CY0720_ALPHA1_BETA1_405.dat')
Cy_r30_old = CY0720_ALPHA1_BETA1_405;
Cy_r30 = reshape(Cy_r30_old,[20,19]);

%% Cn_r30
load('CN0720_ALPHA1_BETA1_503.dat')
Cn_r30_old = CN0720_ALPHA1_BETA1_503;
Cn_r30 = reshape(Cn_r30_old,[20,19]);

%% Cl_r30
load('CL0720_ALPHA1_BETA1_603.dat')
Cl_r30_old = CL0720_ALPHA1_BETA1_603;
Cl_r30 = reshape(Cl_r30_old,[20,19]);

%% Cy_a20
load('CY0620_ALPHA1_BETA1_403.dat')
Cy_a20_old = CY0620_ALPHA1_BETA1_403;
Cy_a20 = reshape(Cy_a20_old,[20,19]);

%% Cy_a20_lef
load('CY0920_ALPHA2_BETA1_404.dat')
Cy_a20_lef_old = CY0920_ALPHA2_BETA1_404;
Cy_a20_lef = reshape(Cy_a20_lef_old,[14,19]);

%% Cn_a20
load('CN0620_ALPHA1_BETA1_504.dat')
Cn_a20_old = CN0620_ALPHA1_BETA1_504;
Cn_a20 = reshape(Cn_a20_old,[20,19]);

%% Cn_a20_lef
load('CN0920_ALPHA2_BETA1_505.dat')
Cn_a20_lef_old = CN0920_ALPHA2_BETA1_505;
Cn_a20_lef = reshape(Cn_a20_lef_old,[14,19]);

%% Cl_a20
load('CL0620_ALPHA1_BETA1_604.dat')
Cl_a20_old = CL0620_ALPHA1_BETA1_604;
Cl_a20 = reshape(Cl_a20_old,[20,19]);

%% Cl_a20_lef
load('CL0920_ALPHA2_BETA1_605.dat')
Cl_a20_lef_old = CL0920_ALPHA2_BETA1_605;
Cl_a20_lef = reshape(Cl_a20_lef_old,[14,19]);

%% delta_Cn_beta
load('CN9999_ALPHA1_brett.dat')
delta_Cn_beta_old = CN9999_ALPHA1_brett;
delta_Cn_beta = reshape(delta_Cn_beta_old,[20,1]);

%% delta_Cl_beta
load('CL9999_ALPHA1_brett.dat')
delta_Cl_beta_old = CL9999_ALPHA1_brett;
delta_Cl_beta = reshape(delta_Cl_beta_old,[20,1]);

%% delta_Cm
load('CM9999_ALPHA1_brett.dat')
delta_Cm_old = CM9999_ALPHA1_brett;
delta_Cm = reshape(delta_Cm_old,[20,1]);

%% eta_el
load('ETA_DH1_brett.dat')
eta_el_old = ETA_DH1_brett;
eta_el = reshape(eta_el_old,[5,1]);

%% write these coefficients to csv files for python interpretation

if write == 1
    
    % coefficient files
    
    writematrix(Cl_dh_neg25,append(csv_path,'hifi_C/Cl_dh_neg25.csv'))
    writematrix(Cl_dh_0,append(csv_path,'hifi_C/Cl_dh_0.csv'))
    writematrix(Cl_dh_pos25,append(csv_path,'hifi_C/Cl_dh_pos25.csv'))
    
    writematrix(Cn_dh_neg25,append(csv_path,'hifi_C/Cn_dh_neg25.csv'))
    writematrix(Cn_dh_0,append(csv_path,'hifi_C/Cn_dh_0.csv'))
    writematrix(Cn_dh_pos25,append(csv_path,'hifi_C/Cn_dh_pos25.csv'))
    
    writematrix(Cm_dh_neg25,append(csv_path,'hifi_C/Cm_dh_neg25.csv'))
    writematrix(Cm_dh_neg10,append(csv_path,'hifi_C/Cm_dh_neg10.csv'))
    writematrix(Cm_dh_0,append(csv_path,'hifi_C/Cm_dh_0.csv'))
    writematrix(Cm_dh_pos10,append(csv_path,'hifi_C/Cm_dh_pos10.csv'))
    writematrix(Cm_dh_pos25,append(csv_path,'hifi_C/Cm_dh_pos25.csv'))
    
    writematrix(Cx_dh_neg25,append(csv_path,'hifi_C/Cx_dh_neg25.csv'))
    writematrix(Cx_dh_neg10,append(csv_path,'hifi_C/Cx_dh_neg10.csv'))
    writematrix(Cx_dh_0,append(csv_path,'hifi_C/Cx_dh_0.csv'))
    writematrix(Cx_dh_pos10,append(csv_path,'hifi_C/Cx_dh_pos10.csv'))
    writematrix(Cx_dh_pos25,append(csv_path,'hifi_C/Cx_dh_pos25.csv'))
    
    writematrix(Cy,append(csv_path,'hifi_C/Cy.csv'))
    
    writematrix(Cz_dh_neg25,append(csv_path,'hifi_C/Cz_dh_neg25.csv'))
    writematrix(Cz_dh_neg10,append(csv_path,'hifi_C/Cz_dh_neg10.csv'))
    writematrix(Cz_dh_0,append(csv_path,'hifi_C/Cz_dh_0.csv'))
    writematrix(Cz_dh_pos10,append(csv_path,'hifi_C/Cz_dh_pos10.csv'))
    writematrix(Cz_dh_pos25,append(csv_path,'hifi_C/Cz_dh_pos25.csv'))
    
    % damping files
    
    writematrix(Cxq,append(csv_path,'hifi_damping/Cxq.csv'))
    writematrix(Cyr,append(csv_path,'hifi_damping/Cyr.csv'))
    writematrix(Cyp,append(csv_path,'hifi_damping/Cyp.csv'))
    writematrix(Czq,append(csv_path,'hifi_damping/Czq.csv'))
    writematrix(Clr,append(csv_path,'hifi_damping/Clr.csv'))
    writematrix(Clp,append(csv_path,'hifi_damping/Clp.csv'))
    writematrix(Cmq,append(csv_path,'hifi_damping/Cmq.csv'))
    writematrix(Cnr,append(csv_path,'hifi_damping/Cnr.csv'))
    writematrix(Cnp,append(csv_path,'hifi_damping/Cnp.csv'))
    
    % c_lef files
    
    writematrix(Cx_lef,append(csv_path,'hifi_C_lef/Cx_lef.csv'))
    writematrix(Cy_lef,append(csv_path,'hifi_C_lef/Cy_lef.csv'))
    writematrix(Cz_lef,append(csv_path,'hifi_C_lef/Cz_lef.csv'))
    writematrix(Cl_lef,append(csv_path,'hifi_C_lef/Cl_lef.csv'))
    writematrix(Cm_lef,append(csv_path,'hifi_C_lef/Cm_lef.csv'))
    writematrix(Cn_lef,append(csv_path,'hifi_C_lef/Cn_lef.csv'))
    
    % damping_lef files
    
    writematrix(delta_Cxq_lef,append(csv_path,'hifi_damping_lef/delta_Cxq_lef.csv'))
    writematrix(delta_Cyr_lef,append(csv_path,'hifi_damping_lef/delta_Cyr_lef.csv'))
    writematrix(delta_Cyp_lef,append(csv_path,'hifi_damping_lef/delta_Cyp_lef.csv'))
    writematrix(delta_Czq_lef,append(csv_path,'hifi_damping_lef/delta_Czq_lef.csv'))
    writematrix(delta_Clr_lef,append(csv_path,'hifi_damping_lef/delta_Clr_lef.csv'))
    writematrix(delta_Clp_lef,append(csv_path,'hifi_damping_lef/delta_Clp_lef.csv'))
    writematrix(delta_Cmq_lef,append(csv_path,'hifi_damping_lef/delta_Cmq_lef.csv'))
    writematrix(delta_Cnr_lef,append(csv_path,'hifi_damping_lef/delta_Cnr_lef.csv'))
    writematrix(delta_Cnp_lef,append(csv_path,'hifi_damping_lef/delta_Cnp_lef.csv'))
    
    % rudder files
    
    writematrix(Cy_r30,append(csv_path,'hifi_rudder/Cy_r30.csv'))
    writematrix(Cn_r30,append(csv_path,'hifi_rudder/Cn_r30.csv'))
    writematrix(Cl_r30,append(csv_path,'hifi_rudder/Cl_r30.csv'))   
    
    % aileron files
    
    writematrix(Cy_a20,append(csv_path,'hifi_ailerons/Cy_a20.csv'))
    writematrix(Cy_a20_lef,append(csv_path,'hifi_ailerons/Cy_a20_lef.csv'))
    writematrix(Cn_a20,append(csv_path,'hifi_ailerons/Cn_a20.csv'))
    writematrix(Cn_a20_lef,append(csv_path,'hifi_ailerons/Cn_a20_lef.csv'))
    writematrix(Cl_a20,append(csv_path,'hifi_ailerons/Cl_a20.csv'))
    writematrix(Cl_a20_lef,append(csv_path,'hifi_ailerons/Cl_a20_lef.csv'))
    
    % other coeff files
    
    writematrix(delta_Cn_beta,append(csv_path,'hifi_other_coeffs/delta_Cn_beta.csv'))
    writematrix(delta_Cl_beta,append(csv_path,'hifi_other_coeffs/delta_Cl_beta.csv'))
    writematrix(delta_Cm,append(csv_path,'hifi_other_coeffs/delta_Cm.csv'))
    writematrix(eta_el,append(csv_path,'hifi_other_coeffs/eta_el.csv'))

    % we are missing delta_Cm_ds, the C file sets it to zero to ignore deep
    % stall
    
end
