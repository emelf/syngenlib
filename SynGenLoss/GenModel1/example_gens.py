import numpy as np
from SynGenLoss.SatModels.SatModel1 import SaturationModel1, SatModelDataClass1
from .DataClasses import GenDataClass1
from .GenModel import GeneratorModel1


# 103 MVA Generator: 
model_data = GenDataClass1() 
model_data.standard_params(Sn_mva=103, V_nom_kV=11.0, cos_phi=0.9, If_nom_A=525.15, Ra=0.00182, Xd=1.059, Xq=0.676, Xp=0.141)
model_data.nominal_losses(V_nom=1.0, Ia_nom=5406.1, If_nom=1065, P_an_kW=187.46, P_sn_kW=89.16, P_fn_kW=173.65, P_brn_kW=2.13, P_exn_kW=15.88, 
                          P_cn_kW=211.92, P_wfn_kW=172.92, P_bn_kW=240.90)

sat_model_data = SatModelDataClass1(bv=1.0, k=1.0308, Cm=0.160, m=7, Ra_70=0.00182, Xd_sat=1.059, Xq_sat=0.676, Xp=0.141)
sat_model = SaturationModel1(sat_model_data)

Gen103MVA = GeneratorModel1(model_data, sat_model)


# 160 MVA Generator: 
model_data = GenDataClass1() 
model_data.standard_params(Sn_mva=160, V_nom_kV=15.0, cos_phi=0.95, If_nom_A=594.0, Ra=0.002322, Xd=0.8, Xq=0.6, Xp=0.18)
model_data.nominal_losses(V_nom=1.0, Ia_nom=5406, If_nom=1047.0, P_an_kW=327.05, P_sn_kW=237.07, P_fn_kW=477.81, P_brn_kW=5.3, P_exn_kW=33.96, 
                          P_cn_kW=539.87, P_wfn_kW=710.47, P_bn_kW=156.17)

sat_model_data = SatModelDataClass1(bv=1.0, k=1.0308, Cm=0.016, m=7, Ra_70=0.002322, Xd_sat=0.8, Xq_sat=0.6, Xp=0.18)
sat_model = SaturationModel1(sat_model_data)

Gen160MVA = GeneratorModel1(model_data, sat_model)