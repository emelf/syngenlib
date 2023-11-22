from dataclasses import dataclass
from typing import Tuple, Optional
from math import sqrt

@dataclass
class GenDataClass: 
    """A dataclass for storing generator model parameters. 

    Sn_mva: Rated apparent power of the generator. [MVA] \n 
    V_nom_kV: Rated nominal voltage of the generator. [kV] \n
    cos_phi_nom: Power factor at nominal operating condition. [.] \n 
    X_d_u: Unsaturated d-axis reactance of the generator. [pu] \n 
    delta_max: Rotor angle at stability limit [rad] \n 
    P_g_min_pu: Minimum active power [pu] \n 
    P_g_max_pu: Maximum active power [pu] \n 
    P_loss_nom_stator_pu: Stator power losses at nominal operating point [pu] \n 
    P_loss_nom_rotor_pu: Rotor power losses at nominal operating point [pu] \n
    P_loss_nom_core_pu: Core power losses at nominal operating point [pu] \n
    P_loss_nom_const_pu: Constant power losses, e.g., friction and windage [pu] \n
    E_q_max: Internal voltage related to the rotor heating limit [pu] \n 
    """

    S_n_mva: float 
    V_nom_kV: float 
    cos_phi_nom: float 
    X_d_u: float 
    delta_max: float 
    P_g_min_pu: float 
    P_g_max_pu: float 
    P_loss_nom_stator_pu: float 
    P_loss_nom_rotor_pu: float 
    P_loss_nom_core_pu: float 
    P_loss_nom_const_pu: float 
    X_q_u: Optional[float] = None
    E_q_max: Optional[float] = None 
    V_g_min: Optional[float] = 0.95 
    V_g_max: Optional[float] = 1.05
    I_g_max: Optional[float] = 1.05

    def __post_init__(self): 
        P_nom = self.cos_phi_nom 
        Q_nom = sqrt(1 - P_nom**2) 
        V_nom = 1.0
        self.E_q_nom = V_nom * sqrt((1.0 + self.X_d_u*Q_nom/V_nom**2)**2 + (self.X_d_u*P_nom/V_nom**2)**2) 
        if self.E_q_max is None: 
            self.E_q_max = self.E_q_nom 
        self.k_If = self.E_q_nom**(-1) 
        if self.X_q_u is None: 
            self.X_q_u = self.X_d_u 


@dataclass
class TrafoDataClass: 
    """A dataclass for storing trafo model parameters.
    S_n_mva: Rated apparent power of the transformer. [MVA] \n 
    V_nom_kV: Rated nominal voltage of the LV side of the transformer. [kV] \n
    V_SCH: Short circuit voltage [pu] \n 
    P_Cu: Copper losses during nominal operation [pu] \n 
    I_E: Excitation current (to shunt) [pu] \n 
    P_Fe: Iron losses (shunt power losses) [pu] \n 
    tap_ratio: Trafo winding ratio between primary and secondary. 
    """
    S_n_mva: float 
    V_nom_kV: float 
    V_SCH: float 
    P_Cu: float 
    I_E: float 
    P_Fe: float
    tap_ratio: Optional[float] = 1.0

    def __post_init__(self): 
        self.R_T = self.P_Cu 
        self.Z_T = self.V_SCH 
        self.X_T = sqrt(self.Z_T**2 - self.R_T**2)
        self.G_Fe = self.P_Fe 
        self.B_mu = sqrt(self.I_E**2 - self.G_Fe**2)
        self.Y_M = self.G_Fe - 1j*self.B_mu

@dataclass
class PowerLossResult: 
    P_g_pu: float 
    Q_g_pu: float 
    I_f_pu: float 
    V_g_pu: float 
    P_loss_stator_pu: float
    P_loss_rotor_pu: float
    P_loss_core_pu: float
    P_loss_const_pu: float

    def __post_init__(self): 
        self.P_loss_tot = self.P_loss_stator_pu + self.P_loss_rotor_pu + self.P_loss_core_pu + self.P_loss_const_pu
        self.eff = self.P_g_pu / (self.P_g_pu + self.P_loss_tot) 

    def get_losses_pu(self) -> float: 
        return self.P_loss_tot

    def get_losses_mw(self, S_base) -> float: 
        """Returns power losses in MW"""
        return S_base*self.P_loss_tot
    
    def get_data_pu(self) -> Tuple[float, float, float, float, float]: 
        """Returns (P_g_pu, Q_g_pu, V_g_pu, I_f_pu, P_loss_pu)"""
        return (self.P_g_pu, self.Q_g_pu, self.V_g_pu, self.I_f_pu, self.P_loss_tot)
    
    def get_data_mw(self, S_base) -> Tuple[float, float, float, float, float]: 
        """Returns (P_g_mw, Q_g_mvar, V_g_pu, I_f_pu, P_loss_mw)"""
        return (self.P_g_pu*S_base, self.Q_g_pu*S_base, self.V_g_pu, self.I_f_pu, self.get_losses_mw(S_base))
    
