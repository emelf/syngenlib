import numpy as np
from SynGenLoss.SatModels.SatModel2 import SaturationModel2, SatModelDataClass2
from .DataClasses import GenDataClass2
from .GenModel import GeneratorModel2

i_0 = 1.0
i_10 = 1.2
i_12 = 1.8

SG10 = i_10/i_0 - 1
SG12 = i_12/(1.2 * i_0) - 1

# 103 MVA Generator: 
model_data = GenDataClass2() 
model_data.standard_params(S_n_mva=103, V_nom_kV=11, cos_phi=0.9, I_f_base=525.15, R_a_nom=0.00182, R_f_nom=None, X_d_u=1.087, 
                           X_q_u=0.676, X_l=0.08)
model_data.nominal_losses(V_nom=1.0, I_fd_pu_nom=1065/525.15, P_sn_kW=276.62, P_rn_kW=175.78, P_cn_kW=211.92, P_const_kW=413.82)

sat_model_data = SatModelDataClass2(X_d_u=1.087, X_q_u=0.676, X_l=0.08, R_a=0.00182, SG10=SG10, SG12=SG12)
sat_model = SaturationModel2(sat_model_data)

Gen103MVA = GeneratorModel2(model_data, sat_model)
