from ..common import GeneratorDataClass 
import numpy as np 

gen_103_mva = GeneratorDataClass(S_n_mva=103, V_nom_kV=11, cos_phi_nom=0.9, 
                                 X_d_u=1.087, X_q_u=0.676, delta_max=30*np.pi/180, 
                                 P_g_min_pu=0.1, P_g_max_pu=1.0, V_g_min=0.95, V_g_max=1.05,
                                 P_loss_nom_stator_pu=(187.46+89.16)/103_000, 
                                 P_loss_nom_rotor_pu=(173.65+15.88+2.13)/103_000, 
                                 P_loss_nom_core_pu=211.92/103_000, 
                                 P_loss_nom_const_pu=(240.9+172.92)/103_000) 

gen_60_mva = GeneratorDataClass(S_n_mva=60, V_nom_kV=9.5, cos_phi_nom=0.86, 
                                X_d_u=1.205889203, X_q_u=0.661362922, R_a=0.002773253, 
                                delta_max=30*np.pi/180, P_g_min_pu=0.1, P_g_max_pu=1.0, 
                                V_g_min=0.95, V_g_max=1.05,
                                P_loss_nom_stator_pu=(187.46+89.16)/103_000, # TODO: Find more repreentative values
                                P_loss_nom_rotor_pu=(173.65+15.88+2.13)/103_000, # TODO: Find more repreentative values 
                                P_loss_nom_core_pu=211.92/103_000, # TODO: Find more repreentative values
                                P_loss_nom_const_pu=(240.9+172.92)/103_000) # TODO: Find more repreentative values

gen_120_mva = GeneratorDataClass(S_n_mva=120, V_nom_kV=14.7, cos_phi_nom=0.86, 
                                 X_d_u=1.217065518, X_q_u=0.633719567, R_a=0.001874355, 
                                 delta_max=30*np.pi/180, P_g_min_pu=0.1, P_g_max_pu=1.0, 
                                 V_g_min=0.95, V_g_max=1.05,
                                 P_loss_nom_stator_pu=(187.46+89.16)/103_000, # TODO: Find more repreentative values
                                 P_loss_nom_rotor_pu=(173.65+15.88+2.13)/103_000, # TODO: Find more repreentative values 
                                 P_loss_nom_core_pu=211.92/103_000, # TODO: Find more repreentative values
                                 P_loss_nom_const_pu=(240.9+172.92)/103_000) # TODO: Find more repreentative values


test_generator = GeneratorDataClass(S_n_mva=100, V_nom_kV=10, cos_phi_nom=0.9, 
                                    X_d_u=1.0, X_q_u=0.5, delta_max=30*np.pi/180, 
                                    P_g_min_pu=0.1, P_g_max_pu=1.0, V_g_min=0.95, V_g_max=1.05,
                                    P_loss_nom_stator_pu=0.01, 
                                    P_loss_nom_rotor_pu=0.01, 
                                    P_loss_nom_core_pu=0.01, 
                                    P_loss_nom_const_pu=0.02) 