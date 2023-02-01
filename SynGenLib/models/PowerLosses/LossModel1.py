from math import sqrt, atan, pi
import cmath as cm
from typing import Sequence, Optional, Tuple
import numpy as np
import numdifftools as nd 

from SynGenLib.models.SatModels.SatModel1 import SaturationModel1 
from ...utils import GenLossRes

class GenDataClass1: 
    """A dataclass for storing generator model parameters. \n 
    Usage: Call the following functions \n 
    standard_params(*params) \n 
    nominal_losses(*params) \n """
    
    def standard_params(self, Sn_mva:float, V_nom_kV: float, cos_phi: float, If_nom_A: float, 
                        Ra: float, Xd: float, Xq: float, Xp: float) -> None: 
        """Sn_mva: Rated apparent power of the generator. [MVA] \n 
        V_nom_kV: Rated nominal voltage of the generator. [kV] \n
        cos_phi: Power factor at nominal operating condition. [.] \n 
        Ia_nom_A: Nominal armature current of the generator. [A] \n
        If_nom_A: Nominal field current at no-load V = 1.0 pu. [A] \n 
        Ra: Armature resistance of the generator. [pu] \n 
        Xd: Direct axis synchronous reactance of the generator. [pu] \n 
        Xq: Quadrature axis synchronous reactance of the generator. [pu] \n
        Xp: Potier reactance of the generator. [pu] \n"""
        self.Sn_mva = Sn_mva
        self.V_nom_kV = V_nom_kV
        self.cos_phi = cos_phi
        self.Ia_nom_A = Sn_mva*1e6/(sqrt(3)*V_nom_kV*1e3)
        self.If_nom_A = If_nom_A
        self.Ra, self.Xd, self.Xq, self.Xp = Ra, Xd, Xq, Xp
        
    def nominal_losses(self, V_nom:float, Ia_nom:float, If_nom:float, P_an_kW:float, P_sn_kW:float, 
                       P_fn_kW:float, P_brn_kW:float, P_exn_kW:float, P_cn_kW:float, P_wfn_kW:float, 
                       P_bn_kW:float) -> None:
        """V_nom: Voltage on generator terminal for given power loss values [pu] \n 
        Ia_nom: Armature current for given power loss values [A] \n 
        If_nom: Field current for given power loss values [A] \n 
        P_an_kW: Copper power loss from the armature [kW] \n 
        P_sn_kW: Stray power loss from the machine [kW] \n 
        P_fn_kW: Field power loss from the machine [kW] \n
        P_brn_kW: Brush power loss from the machine [kW] \n 
        P_exn_kW: Exciter power loss from the machine [kW] \n 
        P_cn_kW: Iron core power loss from the machine [kW] \n 
        P_wfn_kW: Winding and friction power loss from the machine [kW] \n 
        P_bn_kW: Bearing power loss from the machine [kW]"""
        self.V_nom = V_nom
        self.Ia_nom = Ia_nom / self.Ia_nom_A
        self.If_nom = If_nom / self.If_nom_A
        self.P_an = P_an_kW/(self.Sn_mva*1000)
        self.P_sn = P_sn_kW/(self.Sn_mva*1000)
        self.P_fn = P_fn_kW/(self.Sn_mva*1000)
        self.P_brn = P_brn_kW/(self.Sn_mva*1000)
        self.P_exn = P_exn_kW/(self.Sn_mva*1000)
        self.P_cn = P_cn_kW/(self.Sn_mva*1000)
        self.P_wfn = P_wfn_kW/(self.Sn_mva*1000)
        self.P_bn = P_bn_kW/(self.Sn_mva*1000)
        
    def get_standard_params(self) -> Sequence[float]: 
        return self.Sn_mva, self.V_nom_pu, self.cos_phi, self.Ia_nom, self.If_nom, self.Ra, self.Xd, self.Xq, self.Xp
    
    def get_nominal_losses(self) -> Sequence[float]: 
        return self.P_an, self.P_sn, self.P_fn, self.P_brn, self.P_exn, self.P_cn, self.P_wfn, self.P_bn


class GeneratorModel1_If_in: 
    """ Main class for the generator loss model. Requires model data and saturation model to be defined before use. """
    def __init__(self, model_data: GenDataClass1) -> None: 
        self.md = model_data
        self.u_names = ("P_g", "V_g")
        self.y_names = ("Q_g")
    
    def _calc_currents(self, P_pu: float, Q_pu: float, V_t: float) -> float: 
        """Calculates the stator and rotor currents (and load angle) based on given inputs. \n
        returns ia"""
        ia = np.sqrt(P_pu**2 + Q_pu**2)/V_t
        return ia

    def _calc_losses_pu(self, V_t: float, I_a: float, I_f: float) -> Tuple[float, float, float]: 
        """Calculate generator losses based on P, Q, and Vt. \n
        returns an instance of the GenLossRes class. """
        P_loss_stator = (self.md.P_an + self.md.P_sn)*(I_a/self.md.Ia_nom)**2
        P_loss_rotor = (self.md.P_fn + self.md.P_brn)*(I_f/self.md.If_nom)**2
        P_loss_exciter = self.md.P_exn*I_f/self.md.If_nom
        P_loss_core = self.md.P_cn*(V_t/self.md.V_nom)**2
        return P_loss_stator, P_loss_rotor, P_loss_exciter, P_loss_core

    def get_P_losses(self, P_pu: float, Q_pu: float, V_t: float, I_f: float) -> GenLossRes: 
        I_a = self._calc_currents(P_pu, Q_pu, V_t) 
        P_loss_stator, P_loss_rotor, P_loss_exciter, P_loss_core = self._calc_losses_pu(V_t, I_a, I_f)
        gen_loss_res = GenLossRes(P_pu, P_loss_rotor, P_loss_exciter, P_loss_stator, P_loss_core, self.md.P_bn+self.md.P_wfn)
        return gen_loss_res 


class GeneratorModel1(GeneratorModel1_If_in): 
    def __init__(self, model_data: GenDataClass1, sat_model: SaturationModel1): 
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
        P_loss_stator, P_loss_rotor, P_loss_exciter, P_loss_core = self._calc_losses_pu(V_t, I_a, I_f)
        gen_loss_res = GenLossRes(P_pu, P_loss_rotor, P_loss_exciter, P_loss_stator, P_loss_core, self.md.P_bn+self.md.P_wfn)
        return gen_loss_res 