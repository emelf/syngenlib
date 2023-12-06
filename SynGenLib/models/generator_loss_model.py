from typing import Tuple
import numpy as np

from ..common import GeneratorDataClass, GeneratorLossResult

class GeneratorLossModel: 
    def __init__(self, model_data: GeneratorDataClass): 
        self.md = model_data

        P_nom = self.md.cos_phi_nom 
        Q_nom = P_nom * np.tan(np.arccos(P_nom))
        self.E_q_nom = self._calc_E_q(P_nom, Q_nom, 1.0)

    def _calc_E_q(self, P_g, Q_g, V_g):
        E_q_square = V_g**2*((1.0 + self.md.X_d_u*Q_g/(V_g**2))**2 + (self.md.X_d_u*P_g/(V_g**2))**2)
        return np.sqrt(E_q_square)
    
    def _calc_currents(self, P_pu: float, Q_pu: float, V_g: float) -> Tuple[float, float]: 
        """Calculates the stator and rotor currents (and load angle) based on given inputs. \n
        returns (I_a, I_f)"""
        I_a = np.sqrt(P_pu**2 + Q_pu**2)/V_g
        E_q = self._calc_E_q(P_pu, Q_pu, V_g)
        I_f = E_q*self.md.k_If
        return I_a, I_f 

    def _calc_losses_pu(self, V_g: float, I_a: float, I_f: float) -> Tuple[float, float, float]: 
        """Calculate generator losses based on P, Q, and Vt. \n
        returns an instance of the GenLossRes class. """
        P_loss_stator = self.md.P_loss_nom_stator_pu*I_a**2
        P_loss_rotor = self.md.P_loss_nom_rotor_pu*I_f**2
        P_loss_core = self.md.P_loss_nom_core_pu*V_g**2
        return P_loss_stator, P_loss_rotor, P_loss_core

    def get_P_losses(self, P_g_pu: float, Q_g_pu: float, V_g: float) -> GeneratorLossResult: 
        I_a_pu, I_f_pu = self._calc_currents(P_g_pu, Q_g_pu, V_g) 
        P_loss_stator, P_loss_rotor, P_loss_core = self._calc_losses_pu(V_g, I_a_pu, I_f_pu)
        gen_loss_res = GeneratorLossResult(P_g_pu, Q_g_pu, I_f_pu, V_g, P_loss_stator, P_loss_rotor, P_loss_core, self.md.P_loss_nom_const_pu)
        return gen_loss_res 
    
    def get_P_loss_grad_pu(self, P_g_pu: float, Q_g_pu: float, V_g: float) -> Tuple[float, float, float]: 
        """ Returns the power loss gradient w.r.t. P_g, Q_g, V_g respectively. $ 
        """
        dP_s_dP_g = 2*self.md.P_loss_nom_stator_pu*P_g_pu/V_g**2  
        dP_s_dQ_g = 2*self.md.P_loss_nom_stator_pu*Q_g_pu/V_g**2 
        dP_s_dV_g = -(2*self.md.P_loss_nom_stator_pu*(P_g_pu**2 + Q_g_pu**2))/(V_g**3)

        dP_r_dP_g = (2*self.md.P_loss_nom_rotor_pu*P_g_pu*self.md.X_d_u**2 * self.md.k_If**2)/(V_g**2)
        dP_r_dQ_g = 2*self.md.P_loss_nom_rotor_pu*self.md.X_d_u*self.md.k_If**2*(Q_g_pu*self.md.X_d_u/V_g**2 + 1)
        dP_r_dV_g = (2*self.md.P_loss_nom_rotor_pu*self.md.k_If**2*(-P_g_pu**2*self.md.X_d_u**2 - Q_g_pu**2*self.md.X_d_u**2 + V_g**4)) / (V_g**3)

        dP_c_dP_g = 0.0
        dP_c_dQ_g = 0.0
        dP_c_dV_g = 2*self.md.P_loss_nom_core_pu*V_g 

        dP_L_dP_g = dP_s_dP_g + dP_r_dP_g + dP_c_dP_g 
        dP_L_dQ_g = dP_s_dQ_g + dP_r_dQ_g + dP_c_dQ_g
        dP_L_dV_g = dP_s_dV_g + dP_r_dV_g + dP_c_dV_g

        return (dP_L_dP_g, dP_L_dQ_g, dP_L_dV_g) 
