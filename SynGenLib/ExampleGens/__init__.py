from ..DataClasses import GenDataClass 
import numpy as np 

gen_103_mva = GenDataClass(S_n_mva=103, V_nom_kV=11, cos_phi_nom=0.9, 
                           X_d_u=1.087, X_q_u=0.676, delta_max=30*np.pi/180, 
                           P_g_min_pu=0.1, P_g_max_pu=1.0, V_g_min=0.95, V_g_max=1.05,
                           P_loss_nom_stator_pu=(187.46+89.16)/103_000, 
                           P_loss_nom_rotor_pu=(173.65+15.88+2.13)/103_000, 
                           P_loss_nom_core_pu=211.92/103_000, 
                           P_loss_nom_const_pu=(240.9+172.92)/103_000) 

