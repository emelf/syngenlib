import numpy as np
from typing import Optional, Tuple
from SynGenLoss.SatModels.SatModel1 import SaturationModel1 
from .DataClasses import GenDataClass1
from ..utils import GenLossRes

def get_core_losses(Vt: float, gen_data: GenDataClass1) -> float: 
    """Extrapolates the constant losses with respect to the terminal voltage. \n 
    Vt: Terminal generator voltage. \n 
    P_cn: Nominal core losses at Vt = 1.0. \n 
    P_wfn: Windage and friction losses. \n 
    P_bn: Nominal bearings losses. \n 
    Output: Constant losses [pu]"""
    P_core = gen_data.P_cn*(Vt/gen_data.V_nom)**2
    return P_core

def get_rotor_loss(If: float, gen_data: GenDataClass1) -> float: 
    """Extrapolates the rotor copper, brush, and exciter losses from nominal values using the field current If. \n 
    If: Field current of given operating point. \n """
    P_rotor = (gen_data.P_fn + gen_data.P_brn)*(If/gen_data.If_nom)**2
    return P_rotor

def get_stator_loss(Ia: float, gen_data: GenDataClass1) -> float:
    """Extrapolates the stator armature and stray losses based on the armature current. \n 
    Ia: Armature current [pu] \n Ia_nom: Nominal armature current [pu] \n
    P_an: Nominal armature losses [pu]  \n
    P_sn: Nominal stray losses [pu] \n 
    Output: Stator losses [pu]"""
    P_stator = (gen_data.P_an + gen_data.P_sn)*(Ia/gen_data.Ia_nom)**2
    return P_stator

def get_exciter_loss(If:float, gen_data: GenDataClass1) -> float: 
    return gen_data.P_exn*If/gen_data.If_nom
     

class GeneratorModel1_If_in: 
    """ Main class for the generator loss model. Requires model data and saturation model to be defined before use. """
    def __init__(self, model_data: GenDataClass1) -> None: 
        self.md = model_data
    
    def _calc_currents(self, P_pu: float, Q_pu: float, V_t: float) -> float: 
        """Calculates the stator and rotor currents (and load angle) based on given inputs. \n
        returns ia"""
        ia = np.sqrt(P_pu**2 + Q_pu**2)/V_t
        return ia

    def _calc_losses_pu(self, V_t: float, I_a: float, I_f: float) -> Tuple[float, float, float]: 
        """Calculate generator losses based on P, Q, and Vt. \n
        returns an instance of the GenLossRes class. """
        P_loss_stator = get_stator_loss(I_a, self.md)
        P_loss_rotor = get_rotor_loss(I_f, self.md)
        P_loss_exciter = get_exciter_loss(I_f, self.md)
        P_loss_core = get_core_losses(V_t, self.md)
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