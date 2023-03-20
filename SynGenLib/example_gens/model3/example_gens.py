import numpy as np
from ...models.SatModels.SatModel3 import SaturationModel3, SatModelDataClass3

#103 MVA generator
i_0 = 1.0
i_10 = 1.2
i_12 = 1.8

SG10 = i_10/i_0 - 1
SG12 = i_12/(1.2 * i_0) - 1

# 103 MVA Generator: 
gd = {"S_n_mva": 103, 'V_nom_kV': 11, 'cos_phi': 0.9, "I_f_base": 525.15, "R_a_nom": 0.00182, "R_f_nom": None, "X_d_u": 1.059, "X_q_u": 0.676, "X_l": 0.08, 
      "V_nom": 1.0, "I_fd_pu_nom": 1055/525.15, "P_sn_kW": 276.62, "P_rn_kW": 175.78, "P_exn_kW": 15.88, "P_cn_kW": 211.92, "P_const_kW": 413.82, 
      "SG10": SG10, "SG12": SG12, "P_min": 1e-1, "P_max": 1.01, "I_f_max": 1065/525.15}

# model_data_103 = GenDataClass2() 
# model_data_103.standard_params(S_n_mva=gd["S_n_mva"], V_nom_kV=gd["V_nom_kV"], cos_phi=gd["cos_phi"], I_f_base=gd["I_f_base"], R_a_nom=gd["R_a_nom"], 
#                                R_f_nom=gd["R_f_nom"], X_d_u=gd["X_d_u"], X_q_u=gd["X_q_u"], X_l=gd["X_l"])
# model_data_103.nominal_losses(V_nom=gd["V_nom"], I_fd_pu_nom=gd["I_fd_pu_nom"], P_sn_kW=gd["P_sn_kW"], P_rn_kW=gd["P_rn_kW"], P_exn_kW=gd["P_exn_kW"], 
#                               P_cn_kW=gd["P_cn_kW"], P_const_kW=gd["P_const_kW"])

sat_model_data_103 = SatModelDataClass3(X_d_u=gd["X_d_u"], X_q_u=gd["X_q_u"], X_l=gd["X_l"], R_a=gd["R_a_nom"], SG10=gd["SG10"], SG12=gd["SG12"])
sat_model_103 = SaturationModel3(sat_model_data_103)
# GenCD_data = CapabilityDataClass1(P_lims=(gd["P_min"], gd["P_max"]), X_d=gd["X_d_u"], X_q=gd["X_q_u"], R_a=gd["R_a_nom"], I_f_max=gd["I_f_max"])

# Gen103MVA = GeneratorModel2(model_data_103, sat_model_103)
# Gen103MVA_If_in = GeneratorModel2_If_in(model_data_103)
Gen103MVA_sat_model = sat_model_103 
# Gen103MVA_CD = CapabilityDiagram1(GenCD_data, Gen103MVA_sat_model)
# Gen103MVA_CD_approx = CapabilityDiagram1Approx(Gen103MVA_CD, N_points=100)
