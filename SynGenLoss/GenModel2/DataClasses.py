from math import sqrt, atan, pi
import cmath as cm
from typing import Sequence

class GenDataClass2: 
    """A dataclass for storing generator model parameters. \n 
    Usage: Call the following functions \n 
    standard_params(*params) \n 
    nominal_losses(*params) \n """
    
    def standard_params(self, S_n_mva:float, V_nom_kV: float, cos_phi: float, I_f_nom_A: float, 
                        R_a_nom: float, R_f_nom: float, X_d_u: float, X_q_u: float, X_l: float) -> None: 
        """S_n_mva: Rated apparent power of the generator. [MVA] \n 
        V_nom_kV: Rated nominal voltage of the generator. [kV] \n
        cos_phi: Power factor at nominal operating condition. [.] \n 
        I_f_nom_A: Nominal field current at no-load V = 1.0 pu. [A] \n 
        R_a_nom: Armature resistance of the generator. [pu] \n 
        R_f_nom: Rotor field resistance [pu] \n
        X_d_u: Unsaturated direct axis synchronous reactance of the generator. [pu] \n 
        X_q_u: Unsaturated quadrature axis synchronous reactance of the generator. [pu] \n 
        X_l: Leakage reactance of the generator [pu] """
        self.Sn_mva = S_n_mva
        self.V_nom_kV = V_nom_kV
        self.cos_phi = cos_phi
        self.Ia_nom_A = S_n_mva*1e6/(sqrt(3)*V_nom_kV*1e3)
        self.If_nom_A = I_f_nom_A
        self.R_a_nom, self.R_f_nom, self.X_d_u, self.X_q_u, self.X_l = R_a_nom, R_f_nom, X_d_u, X_q_u, X_l
        self.X_ad_u = self.X_d_u - self.X_l 
        self.X_aq_u = self.X_q_u - self.X_l

    def nominal_losses(self, I_fd_pu_nom: float, P_sn_kW: float, P_rn_kW: float,
                       P_cn_kW: float, P_const_kW: float): 
        """
        P_sn_kW: Power losses [kW] in stator (+ stray) at nominal oeprating point. \n 
        P_rn_kW: Power losses [kW] in rotor (field + brushes) at nominal operating point. \n 
        P_cn_kW: Power losses [kW] in the core at nominal voltage \n 
        P_const_kW: Power losses [kW] that is constant (bearing + friction and windage) """
        self.P_sn_kW = P_sn_kW 
        self.P_rn_kW = P_rn_kW
        self.P_cn_kW = P_cn_kW
        self.P_const_kW = P_const_kW 

        self.R_st = (self.P_sn_kW/self.Sn_mva/1000) # Assumes I_a = 1.0 pu 
        self.R_rt = (self.P_rn_kW/self.Sn_mva/1000)/I_fd_pu_nom**2 

    def get_standard_params(self) -> Sequence[float]: 
        return self.Sn_mva, self.V_nom_pu, self.cos_phi, self.Ia_nom, self.If_nom, self.Ra, self.Xd, self.Xq, self.Xp
    
    def get_nominal_losses(self) -> Sequence[float]: 
        return self.P_an, self.P_sn, self.P_fn, self.P_brn, self.P_exn, self.P_cn, self.P_wfn, self.P_bn

