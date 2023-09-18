from dataclasses import dataclass
from typing import Tuple

@dataclass
class GenDataClass: 
    """A dataclass for storing generator model parameters. 

    Sn_mva: Rated apparent power of the generator. [MVA] \n 
    V_nom_kV: Rated nominal voltage of the generator. [kV] \n
    cos_phi_nom: Power factor at nominal operating condition. [.] \n 
    X_d: Direct axis synchronous reactance of the generator. [pu] \n 
    E_q_max: Internal voltage related to the rotor heating limit [pu] \n 
    delta_max: Rotor angle at stability limit [rad] \n 
    P_g_min_pu: Minimum active power [pu] \n 
    P_g_max_pu: Maximum active power [pu] \n 
    P_loss_nom_stator_pu: Stator power losses at nominal operating point [pu] \n 
    P_loss_nom_rotor_pu: Rotor power losses at nominal operating point [pu] \n
    P_loss_nom_core_pu: Core power losses at nominal operating point [pu] \n
    P_loss_nom_const_pu: Constant power losses, e.g., friction and windage [pu] \n
    """

    S_n_mva: float 
    V_nom_kV: float 
    cos_phi_nom: float 
    X_d: float 
    E_q_max: float 
    delta_max: float 
    P_g_min_pu: float 
    P_g_max_pu: float 
    P_loss_nom_stator_pu: float 
    P_loss_nom_rotor_pu: float 
    P_loss_nom_core_pu: float 
    P_loss_nom_const_pu: float 


class GenLossRes: 
    def __init__(self, P_g_pu: float, Q_g_pu: float, I_f_pu: float, V_g_pu: float, 
                 P_loss_stator_pu: float, P_loss_rotor_pu: float, P_loss_core_pu: float, P_loss_const_pu: float): 
        self.P_g = P_g_pu
        self.Q_g = Q_g_pu
        self.I_f = I_f_pu
        self.V_g = V_g_pu

        self.P_stator_pu = P_loss_stator_pu 
        self.P_rotor_pu = P_loss_rotor_pu 
        self.P_core_pu = P_loss_core_pu 
        self.P_const_pu = P_loss_const_pu 

        self.P_loss_tot = P_loss_stator_pu + P_loss_rotor_pu + P_loss_core_pu + P_loss_const_pu
        self.eff = P_g_pu / (P_g_pu + self.P_loss_tot) 

    def get_losses_pu(self) -> float: 
        return self.P_loss_tot

    def get_losses_mw(self, S_base) -> float: 
        """Returns power losses in MW"""
        return S_base*self.P_loss_tot
    
    def get_data_pu(self) -> Tuple[float, float, float, float, float]: 
        """Returns (P_g_pu, Q_g_pu, V_g_pu, I_f_pu, P_loss_pu)"""
        return (self.P_g, self.Q_g, self.V_g, self.I_f, self.P_loss_tot)
    
    def get_data(self, S_base) -> Tuple[float, float, float, float, float]: 
        """Returns (P_g_mw, Q_g_mvar, V_g_pu, I_f_pu, P_loss_mw)"""
        return (self.P_g*S_base, self.Q_g*S_base, self.V_g, self.I_f, self.get_losses_mw(S_base))
    
