from dataclasses import dataclass
from typing import Tuple, Optional, Sequence
from math import sqrt
import numpy as np

@dataclass
class GeneratorDataClass: 
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
    E_q_min: Minimum internal voltage supported by the excitation system [pu] \n 
    V_g_min: Minimum allowed generator terminal voltage [pu] \n
    V_g_max: Maximum allowed generator terminal voltage [pu] \n
    I_g_max: Maximum steady-state stator current without overloading the generator [pu] \n 
    R_a: Resistance at nominal operating point (e.g., 75 deg. C) [pu]  
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
    E_q_min: Optional[float] = 0.1
    V_g_min: Optional[float] = 0.95 
    V_g_max: Optional[float] = 1.05
    I_g_max: Optional[float] = 1.0
    R_a: Optional[float] = 0.0

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
class TransformerDataClass: 
    """A dataclass for storing trafo model parameters.
    S_n_mva: Rated apparent power of the transformer. [MVA] \n 
    V_nom_kV: Rated nominal voltage of the LV side of the transformer. [kV] \n
    V_SCH: Short circuit voltage [pu] \n 
    P_Cu: Copper losses during nominal operation [pu] \n 
    I_E: Excitation current (to shunt) [pu] \n 
    P_Fe: Iron losses (shunt power losses) [pu] \n 
    tap_ratio: Trafo winding ratio between primary and secondary. \n 
    Z_lv_ratio: How much of the series impedance is located at the LV side (generator side) 
    """
    S_n_mva: float 
    V_nom_kV: float 
    V_SCH: float 
    P_Cu: float 
    I_E: float 
    P_Fe: float
    tap_ratio: Optional[float] = 1.0
    Z_lv_ratio: Optional[float] = 0.5

    def __post_init__(self): 
        self.R_T = self.P_Cu 
        self.X_T = sqrt(self.V_SCH**2 - self.R_T**2)
        self.Z_T = self.R_T + 1j*self.X_T
        self.G_Fe = self.P_Fe 
        self.B_mu = sqrt(self.I_E**2 - self.G_Fe**2)
        self.Y_M = self.G_Fe - 1j*self.B_mu

        self._calc_pi_parameters() 
        self._calc_power_matrix()
        
    def _calc_pi_parameters(self): 
        # Step 1: Get all parameters for the model 
        Y_lv = 1/(self.Z_T * self.Z_lv_ratio)
        Y_hv = 1/(self.Z_T * (1.0-self.Z_lv_ratio)) 

        # Account for the tap changer 
        Y_hv_12 = Y_hv*self.tap_ratio 
        Y_hv_1 = Y_hv*(1-self.tap_ratio)
        Y_hv_2 = Y_hv*self.tap_ratio*(self.tap_ratio - 1)
        Y_lv_12 = Y_lv*self.tap_ratio**(-2)

        # Collect and transform from star to delta 
        Y_1_star = Y_hv_12
        Y_2_star = Y_lv_12 
        Y_3_star = Y_hv_2 + self.Y_M #/2

        Y_num = Y_1_star + Y_2_star + Y_3_star 
        Y_12 = Y_1_star * Y_2_star / Y_num 
        Y_23 = Y_2_star * Y_3_star / Y_num 
        Y_31 = Y_3_star * Y_1_star / Y_num 
        self.Y_hv = Y_31 + Y_hv_1
        self.Y_lv = Y_23 #+ Y_m/2 * self.b1**2
        self.Z_12 = Y_12**-1

    def _calc_power_matrix(self): 
        """Used to calculate voltage and current at the sending side, given recieving values of V_r and I_r. 
        [V_s; I_s] = [A, B; C, D] [V_r; I_r], which is noted: 
        [V_s; I_s] = ABCD_mat [V_r; I_r]"""
        self.A = self.Y_hv*self.Z_12 + 1 
        self.B = self.Z_12 
        self.C = self.Y_lv*self.Y_hv*self.Z_12 + self.Y_lv + self.Y_hv 
        self.D = self.Y_lv*self.Z_12 + 1 
        self.ABCD_mat = np.array([[self.A, self.B], [self.C, self.D]])

    def change_base(self, S_new_mva: float, V_new_kV: float, inplace: Optional[bool]=False):
        """Most common use-case is when the generator and transformer is different ratings. It is then required to 
        convert the transformer to the same base units as the generator for correct calculations. 
        S_new_mva: New base power [MVA] 
        V_new_kV: New base voltage [kV] 
        inplace: True if the dataclass is changed in-place. If false, returns a new transformer object""" 
        Z_b_old = self.V_nom_kV**2/self.S_n_mva 
        Z_b_new = V_new_kV**2/S_new_mva 
        Z_change = Z_b_old/Z_b_new 
        Y_change = 1/Z_change 
        V_SCH_new = self.V_SCH*Z_change  
        P_Cu_new = self.P_Cu*Z_change 
        I_E_new = self.I_E*Y_change 
        P_Fe_new = self.P_Fe*Y_change 
        if inplace: 
            self.V_SCH = V_SCH_new 
            self.P_Cu = P_Cu_new 
            self.I_E = I_E_new 
            self.P_Fe = P_Fe_new
            self.__post_init__() 
        else: 
            return TransformerDataClass(S_new_mva, V_new_kV, V_SCH_new, P_Cu_new, 
                                        I_E_new, P_Fe_new, self.tap_ratio) 
        
    def change_tap_ratio(self, new_tap_ratio): 
        self.tap_ratio = new_tap_ratio 
        self.__post_init__()


@dataclass
class GeneratorLossResult: 
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
    

@dataclass
class TrafoLossResult: # TODO
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


@dataclass
class CapabilityResult: 
    """ limiter_Q_lim: {
    0 => stator limit, 
    1 => rotor limit, 
    2 => stability limit, 
    3 => voltage limit
    }"""
    P_vals: Sequence[float]
    Q_min_tot: Sequence[float]
    Q_max_tot: Sequence[float]
    Q_stator_min: Sequence[float]
    Q_stator_max: Sequence[float]
    Q_rotor_min: Sequence[float]
    Q_rotor_max: Sequence[float]
    Q_stab_min: Sequence [float]
    Q_v_min: Sequence[float]
    Q_v_max: Sequence[float]
    valid_stator_current: Sequence[bool]
    valid_rotor_current: Sequence[bool]
    valid_active_power: Sequence[bool]
    valid_voltage_levels: Sequence[bool]
    limiter_Q_min: Sequence[int]
    limiter_Q_max: Sequence[int]
    
    def get_PQ_plot(self): 
        """Returns (P_plot, Q_plot, valid_stator_current, valid_active_power, valid_voltage_levels)"""
        Q_cap_plot = np.concatenate([self.Q_min_tot, self.Q_max_tot[::-1]])
        P_cap_plot = np.concatenate([self.P_vals, self.P_vals[::-1]])
        valid_stator_current_plot = np.concatenate([self.valid_stator_current, self.valid_stator_current[::-1]])
        valid_active_power_plot = np.concatenate([self.valid_active_power, self.valid_active_power[::-1]])
        valid_voltage_levels_plot = np.concatenate([self.valid_voltage_levels, self.valid_voltage_levels[::-1]])
        return P_cap_plot, Q_cap_plot, valid_stator_current_plot, valid_active_power_plot, valid_voltage_levels_plot
