from typing import Tuple
import numpy as np

from model_dataclass import GenDataClass, GenLossRes    

class GeneratorLossModel: 
    def __init__(self, model_data: GenDataClass) -> None: 
        self.md = model_data

        P_nom = self.md.cos_phi_nom 
        Q_nom = P_nom * np.tan(np.arccos(P_nom))
        self.E_q_nom = self._calc_E_q(P_nom, Q_nom)

    def _calc_E_q(self, P_g, Q_g, V_g): # TODO
        E_q = np.sqrt((1.0 + self.md.X_d*Q_g)**2 + (self.md.X_d*P_g)**2)
        return E_q
    
    def _calc_currents(self, P_pu: float, Q_pu: float, V_t: float) -> float: 
        """Calculates the stator and rotor currents (and load angle) based on given inputs. \n
        returns ia"""
        I_a = np.sqrt(P_pu**2 + Q_pu**2)/V_t
        E_q = self._calc_E_q(P_pu, Q_pu)
        I_f = E_q/self.E_q_nom 
        return I_a, I_f 

    def _calc_losses_pu(self, V_t: float, I_a: float, I_f: float) -> Tuple[float, float, float]: 
        """Calculate generator losses based on P, Q, and Vt. \n
        returns an instance of the GenLossRes class. """
        P_loss_stator = self.md.P_loss_nom_stator_pu*I_a**2
        P_loss_rotor = self.md.P_loss_nom_rotor_pu*I_f**2
        P_loss_core = self.md.P_loss_nom_core_pu*V_t**2
        return P_loss_stator, P_loss_rotor, P_loss_core

    def get_P_losses(self, P_g_pu: float, Q_g_pu: float, V_g: float) -> GenLossRes: 
        I_a_pu, I_f_pu = self._calc_currents(P_g_pu, Q_g_pu, V_g) 
        P_loss_stator, P_loss_rotor, P_loss_core = self._calc_losses_pu(V_g, I_a_pu, I_f_pu)
        gen_loss_res = GenLossRes(P_g_pu, Q_g_pu, I_f_pu, V_g, P_loss_stator, P_loss_rotor, P_loss_core, self.md.P_loss_nom_const_pu)
        return gen_loss_res 


