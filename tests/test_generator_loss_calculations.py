import unittest
from syngenlib.data import GeneratorDataclass, GeneratorLossDataclass
from syngenlib.models import GeneratorCalculationModel, LinearSaturationModel
from syngenlib.data import GeneratorOperatingPoint, GeneratorBranchResults

gen_data = GeneratorDataclass(S_n_mva=100.0, V_nom_kV=11.0, cos_phi=0.9, 
                              X_d_u=0.9, X_q_u=0.8, R_a=0.01) 
pl_data = GeneratorLossDataclass(P_loss_nom_stator_pu=0.01, 
                                 P_loss_nom_rotor_pu=0.01, 
                                 P_loss_nom_core_pu=0.03, 
                                 P_loss_nom_const_pu=0.03)
gen_model = GeneratorCalculationModel(gen_data, power_loss_data=pl_data)  

op1 = GeneratorOperatingPoint(P_mw=10, Q_mvar=0.0, V_pu=1.0)
op2 = GeneratorOperatingPoint(P_mw=90, Q_mvar=-50.0, V_pu=1.0)
op3 = GeneratorOperatingPoint(P_mw=90, Q_mvar=50.0, V_pu=1.0)
op4 = GeneratorOperatingPoint(P_mw=90, Q_mvar=-50.0, V_pu=0.9)
op5 = GeneratorOperatingPoint(P_mw=90, Q_mvar=50.0, V_pu=1.1)

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

def get_E_q_from_op(op: GeneratorOperatingPoint, gen_model: GeneratorCalculationModel): 
    P_g_pu = op.P_mw / gen_model.gen_data.S_n_mva
    Q_g_pu = op.Q_mvar / gen_model.gen_data.S_n_mva
    V_g_pu = op.V_pu

    E_q = gen_model._get_E_q(P_g_pu, Q_g_pu, V_g_pu)
    br = GeneratorBranchResults(op.P_mw, op.Q_mvar, op.P_mw, op.Q_mvar, op.V_pu, op.V_pu)
    I_f = gen_model.saturation_model.get_field_current(br, E_q) 
    return E_q, I_f 

print("E_q, I_f for op1: ", get_E_q_from_op(op1, gen_model))
print("E_q, I_f for op2: ", get_E_q_from_op(op2, gen_model))
print("E_q, I_f for op3: ", get_E_q_from_op(op3, gen_model))
print("E_q, I_f for op4: ", get_E_q_from_op(op4, gen_model))
print("E_q, I_f for op5: ", get_E_q_from_op(op5, gen_model))
