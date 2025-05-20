import unittest
from syngenlib.data import GeneratorDataclass, GeneratorLossDataclass
from syngenlib.models import BranchCalculationModel, LinearSaturationModel
from syngenlib.data import GeneratorOperatingPoint, GeneratorBranchResults

gen_data = GeneratorDataclass(S_n_mva=100.0, V_nom_kV=11.0, cos_phi=0.9, 
                              X_d_u=0.9, X_q_u=0.8, R_a=0.01) 
pl_data = GeneratorLossDataclass(P_loss_nom_stator_pu=0.01, 
                                 P_loss_nom_rotor_pu=0.01, 
                                 P_loss_nom_core_pu=0.03, 
                                 P_loss_nom_const_pu=0.03)
gen_model = BranchCalculationModel(gen_data, power_loss_data=pl_data)  

op1 = GeneratorOperatingPoint(P_mw=10, Q_mvar=0.0, V_kv=1.0 * gen_data.V_nom_kV)
op2 = GeneratorOperatingPoint(P_mw=90, Q_mvar=-50.0, V_kv=1.0 * gen_data.V_nom_kV)
op3 = GeneratorOperatingPoint(P_mw=90, Q_mvar=50.0, V_kv=1.0 * gen_data.V_nom_kV)
op4 = GeneratorOperatingPoint(P_mw=90, Q_mvar=-50.0, V_kv=0.9 * gen_data.V_nom_kV)
op5 = GeneratorOperatingPoint(P_mw=90, Q_mvar=50.0, V_kv=1.1 * gen_data.V_nom_kV)

res_1 = gen_model.get_branch_losses(op1)
res_2 = gen_model.get_branch_losses(op2)
res_3 = gen_model.get_branch_losses(op3)
res_4 = gen_model.get_branch_losses(op4)
res_5 = gen_model.get_branch_losses(op5)

print(res_1)
print(res_2)
print(res_3)
print(res_4)
print(res_5)
