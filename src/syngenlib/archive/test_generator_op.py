from syngenlib.data import GeneratorDataclass, TransformerDataclass, CapabilityModelDataclass, GeneratorLossDataclass
from syngenlib.models import GeneratorCalculationModel, LinearSaturationModel 
from syngenlib.data import GeneratorOperatingPoint, BranchOperatingPoint, PlantOperatingPoint, GeneratorBranchResults

gen_data = GeneratorDataclass(S_n_mva=30.0, V_nom_kV=9.5, cos_phi=0.9, X_d_u=1.1, X_q_u=0.8, R_a=0.003)
sat_model = LinearSaturationModel(gen_data)
capability_data = CapabilityModelDataclass.default_limits(gen_data)

"""
List of transformer properties: 
T1: 30 MVA, 9.5 - 23 kV, tap 1.0, no magnetizing branch, 
T2: 30 MVA, 9.5 - 23 kV, tap 1.0, with magnetizing branch
T3: 30 MVA, 9.5 - 23 kV, tap 0.9, with magnetizing branch,
T4: 30 MVA, 9.5 - 23 kV, tap 1.1, with magnetizing branch
T5: 50 MVA, 9.5 - 23 kV, tap 1.0, with magnetizing branch
T6: 50 MVA, 9.5 - 23 kV, tap 1.1, with magnetizing branch
T7: 30 MVA, 10.0 - 23 kV, tap 1.0, with magnetizing branch
T8: 30 MVA, 9.5 - 25 kV, tap 1.0, with magnetizing branch
T9: 30 MVA, 9.5 - 25 kV, tap 0.9, with magnetizing branch
"""
T1 = TransformerDataclass(30.0, 9.5, 23.0, 0.1, 0.01, 0.0, 0.0, 1.0, 0.5)
T2 = TransformerDataclass(30.0, 9.5, 23.0, 0.1, 0.01, 0.003, 0.001, 1.0, 0.5)
T3 = TransformerDataclass(30.0, 9.5, 23.0, 0.1, 0.01, 0.003, 0.001, 0.9, 0.5)
T4 = TransformerDataclass(30.0, 9.5, 23.0, 0.1, 0.01, 0.003, 0.001, 1.1, 0.5)
T5 = TransformerDataclass(50.0, 9.5, 23.0, 0.1, 0.01, 0.003, 0.001, 1.0, 0.5)
T6 = TransformerDataclass(50.0, 9.5, 23.0, 0.1, 0.01, 0.003, 0.001, 1.1, 0.5)
T7 = TransformerDataclass(30.0, 10.0, 23.0, 0.1, 0.01, 0.003, 0.001, 1.0, 0.5)
T8 = TransformerDataclass(30.0, 9.5, 25.0, 0.1, 0.01, 0.003, 0.001, 1.0, 0.5)
T9 = TransformerDataclass(30.0, 9.5, 25.0, 0.1, 0.01, 0.003, 0.001, 0.9, 0.5)

op1 = GeneratorOperatingPoint(P_mw=25.0, Q_mvar=0.0, V_pu=1.004882855) 
op2 = GeneratorOperatingPoint(P_mw=25.0, Q_mvar=0.0, V_pu=1.004738249) 
op3 = GeneratorOperatingPoint(P_mw=25.0, Q_mvar=0.0, V_pu=1.115930630) 
op4 = GeneratorOperatingPoint(P_mw=25.0, Q_mvar=0.0, V_pu=0.913539313) 
op5 = GeneratorOperatingPoint(P_mw=25.0, Q_mvar=0.0, V_pu=1.003607484) 
op6 = GeneratorOperatingPoint(P_mw=25.0, Q_mvar=0.0, V_pu=0.912801473) 
op7 = GeneratorOperatingPoint(P_mw=25.0, Q_mvar=0.0, V_pu=1.057619209) 
op8 = GeneratorOperatingPoint(P_mw=25.0, Q_mvar=0.0, V_pu=0.924498765) 
op9 = GeneratorOperatingPoint(P_mw=25.0, Q_mvar=0.0, V_pu=1.026995313) 

sol_1 = GeneratorBranchResults(25.0, 0.0, 24.793686423, -2.052793948, 1.004882855, 1.0)
sol_2 = GeneratorBranchResults(25.0, 0.0, 24.763759045, -2.136307742, 1.004738249, 1.0)
sol_3 = GeneratorBranchResults(25.0, 0.0, 24.795781611, -1.76737924, 1.115930630, 1.0)
sol_4 = GeneratorBranchResults(25.0, 0.0, 24.725730924, -2.551981580, 0.913539313, 1.0)
sol_5 = GeneratorBranchResults(25.0, 0.0, 24.825984243, -1.374348229, 1.003607484, 1.0)
sol_6 = GeneratorBranchResults(25.0, 0.0, 24.808756171, -1.607683499, 0.912801473, 1.0)
sol_7 = GeneratorBranchResults(25.0, 0.0, 24.763759044, -2.136307742, 1.057619209, 1.0)
sol_8 = GeneratorBranchResults(25.0, 0.0, 24.731010813, -2.495154196, 0.924498765, 1.0)
sol_9 = GeneratorBranchResults(25.0, 0.0, 24.771259162, -2.052037273, 1.026995313, 1.0)

gm_1 = GeneratorCalculationModel(gen_data, T1, None, sat_model, capability_data)
gm_2 = GeneratorCalculationModel(gen_data, T2, None, sat_model, capability_data)
gm_3 = GeneratorCalculationModel(gen_data, T3, None, sat_model, capability_data)
gm_4 = GeneratorCalculationModel(gen_data, T4, None, sat_model, capability_data)
gm_5 = GeneratorCalculationModel(gen_data, T5, None, sat_model, capability_data)
gm_6 = GeneratorCalculationModel(gen_data, T6, None, sat_model, capability_data)
gm_7 = GeneratorCalculationModel(gen_data, T7, None, sat_model, capability_data)
gm_8 = GeneratorCalculationModel(gen_data, T8, None, sat_model, capability_data)
gm_9 = GeneratorCalculationModel(gen_data, T9, None, sat_model, capability_data)

res_1 = gm_1._calculate_branch_results_from_gen_op(op1)
res_2 = gm_2._calculate_branch_results_from_gen_op(op2)
res_3 = gm_3._calculate_branch_results_from_gen_op(op3)
res_4 = gm_4._calculate_branch_results_from_gen_op(op4)
res_5 = gm_5._calculate_branch_results_from_gen_op(op5)
res_6 = gm_6._calculate_branch_results_from_gen_op(op6)
res_7 = gm_7._calculate_branch_results_from_gen_op(op7)
res_8 = gm_8._calculate_branch_results_from_gen_op(op8)
res_9 = gm_9._calculate_branch_results_from_gen_op(op9)

def print_results(res: GeneratorBranchResults):
    print(f"V_n_calc: {res.V_grid_pu:.6f} V. Q_n = {res.Q_branch_mvar:.6f} Mvar") 

if __name__ == "__main__":
    print("Comparing Powerfactory vs GenOP:")
    print_results(res_1)
    print_results(res_2)
    print_results(res_3)
    print_results(res_4)
    print_results(res_5)
    print_results(res_6)
    print_results(res_7)
    print_results(res_8)
    print_results(res_9)

