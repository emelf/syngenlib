from math import sqrt, atan, pi
import cmath as cm
from typing import Sequence, Optional, Tuple
import numpy as np

from SynGenLib.models.SatModels.SatModel2 import SaturationModel2 
from ...utils import GenLossRes

class GenDataClass2: 
    """A dataclass for storing generator model parameters. \n 
    Usage: Call the following functions \n 
    standard_params(*params) \n 
    nominal_losses(*params) \n """
    
    def standard_params(self, S_n_mva:float, V_nom_kV: float, cos_phi: float, I_f_base: float, R_a_nom: float, R_f_nom: float, 
                        X_d_u: float, X_q_u: float, X_l: float) -> None: 
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
        self.I_a_nom_A = S_n_mva*1e6/(sqrt(3)*V_nom_kV*1e3)
        self.I_f_base = I_f_base
        self.R_a_nom, self.R_f_nom, self.X_d_u, self.X_q_u, self.X_l = R_a_nom, R_f_nom, X_d_u, X_q_u, X_l
        self.X_ad_u = self.X_d_u - self.X_l 
        self.X_aq_u = self.X_q_u - self.X_l

    def nominal_losses(self, V_nom: float, I_fd_pu_nom: float, P_sn_kW: float, P_rn_kW: float, P_exn_kW: float,
                       P_cn_kW: float, P_const_kW: float): 
        """
        V_nom: The voltage at which the losses were estimated. 
        P_sn_kW: Power losses [kW] in stator (+ stray) at nominal oeprating point. \n 
        P_rn_kW: Power losses [kW] in rotor (field + brushes) at nominal operating point. \n 
        P_exn_kW: Power losses [kW] in the static exciter (set to 0 if unknown) \n
        P_cn_kW: Power losses [kW] in the core at nominal voltage \n 
        P_const_kW: Power losses [kW] that is constant (bearing + friction and windage) """
        self.V_nom = V_nom
        self.I_fd_pu_nom = I_fd_pu_nom
        self.P_sn = P_sn_kW / self.Sn_mva / 1000 # Convert to pu 
        self.P_rn = P_rn_kW / self.Sn_mva / 1000
        self.P_exn = P_exn_kW / self.Sn_mva / 1000
        self.P_cn = P_cn_kW / self.Sn_mva / 1000
        self.P_const = P_const_kW  / self.Sn_mva / 1000
        
        self.R_st = self.P_sn # Assumes I_a = 1.0 pu 
        self.R_rt = self.P_rn/I_fd_pu_nom**2 
      

class GeneratorModel2_If_in: 
    """ Main class for the generator loss model. Requires the model data. For loss calculations, field current need to be supplied. """
    def __init__(self, model_data: GenDataClass2) -> None: 
        self.md = model_data
    
    def _calc_currents(self, P_pu: float, Q_pu: float, V_t: float) -> float: 
        """Calculates the stator current. \n
        returns I_a"""
        I_a = np.sqrt(P_pu**2 + Q_pu**2)/V_t
        return I_a

    def _calc_losses_pu(self, V_t: float, I_a: float, I_f: float) -> Tuple[float, float, float, float]: 
        """Calculate generator losses based on V_t, I_a, and I_f. \n
        returns (P_L_stator, P_L_rotor, P_L_ex, P_L_core) """
        P_loss_stator = self.md.R_st*I_a**2
        P_loss_rotor = self.md.R_rt * I_f**2
        P_loss_core = self.md.P_cn*(V_t/self.md.V_nom)**2
        P_loss_ex = I_f/self.md.I_fd_pu_nom * self.md.P_exn
        return P_loss_stator, P_loss_rotor, P_loss_ex, P_loss_core

    def get_P_losses(self, P_pu: float, Q_pu: float, V_t: float, I_f: float) -> GenLossRes: 
        """Based on generator operating point and the excitation current, calculates generator power losses. \n
        returns an instance of the GenLossRes class"""
        I_a = self._calc_currents(P_pu, Q_pu, V_t) 
        P_loss_stator, P_loss_rotor, P_loss_ex, P_loss_core = self._calc_losses_pu(V_t, I_a, I_f)
        gen_loss_res = GenLossRes(P_pu, P_loss_rotor, P_loss_ex, P_loss_stator, P_loss_core, self.md.P_const)
        return gen_loss_res 


class GeneratorModel2(GeneratorModel2_If_in): 
    """If the field current can be calculated from P, Q, V, through the sat_model, then use this class. \n 
    Field current should no longer be supplied to the self.get_P_losses function. """
    def __init__(self, model_data: GenDataClass2, sat_model: SaturationModel2): 
        self.sm = sat_model
        super().__init__(model_data)
    
    def _calc_phi(self, P_el:float, Q_el:float) -> float:
        if P_el == 0 and Q_el == 0: 
            return 0
        elif P_el == 0 and not Q_el == 0: 
            return np.pi/2 * np.sign(Q_el)
        else: 
            return np.arctan(Q_el/P_el) 

    def _calc_currents(self, P_pu: float, Q_pu: float, V_t: float) -> Tuple[float, float, float]: 
        """Calculates the stator and rotor currents (and load angle) based on given inputs. \n
        returns ia [pu], ifd [pu], delta [rad]"""
        ia = np.sqrt(P_pu**2 + Q_pu**2)/V_t
        if hasattr(ia, "__len__"): #Checks if ia is a list/array or a scalar 
            phi = np.array([self._calc_phi(P_el, Q_el) for P_el, Q_el in zip(P_pu, Q_pu)])
        else: 
            phi = self._calc_phi(P_pu, Q_pu)
        ifd, delta, _ = self.sm.calc_ifd(V_t, ia, phi)
        return ia, ifd, delta

    def get_P_losses(self, P_pu: float, Q_pu: float, V_t: float) -> GenLossRes:
        I_a, I_f, _ = self._calc_currents(P_pu, Q_pu, V_t)
        P_loss_stator, P_loss_rotor, P_loss_ex, P_loss_core = self._calc_losses_pu(V_t, I_a, I_f)
        gen_loss_res = GenLossRes(P_pu, P_loss_rotor, P_loss_ex, P_loss_stator, P_loss_core, self.md.P_const)
        return gen_loss_res 