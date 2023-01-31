from math import sqrt, atan, pi
import cmath as cm
from typing import Sequence, Optional, Tuple
import numpy as np
import sympy as sp
from sympy import lambdify, symbols, diff
from dataclasses import dataclass

from SynGenLib.models.SatModels.SatModel2 import SaturationModel2 
from ...utils import GenLossRes

@dataclass
class GenDataClass2: 
    """
    A dataclass for storing generator model parameters. 
    ----------
    S_n_mva: Rated apparent power of the generator. [MVA] 
    V_nom_kV: Rated nominal voltage of the generator. [kV] 
    cos_phi: Power factor at nominal operating condition. [.] 
    I_f_nom_A: Nominal field current at no-load V = 1.0 pu. [A] 
    R_a_nom: Armature resistance of the generator. [pu] 
    R_f_nom: Rotor field resistance [pu] 
    X_d_u: Unsaturated direct axis synchronous reactance of the generator. [pu]
    X_q_u: Unsaturated quadrature axis synchronous reactance of the generator. [pu]
    X_l: Leakage reactance of the generator [pu]
    V_nom: The voltage at which the losses were estimated.
    P_sn_kW: Power losses [kW] in stator (+ stray) at nominal oeprating point.
    P_rn_kW: Power losses [kW] in rotor (field + brushes) at nominal operating point.
    P_exn_kW: Power losses [kW] in the static exciter (set to 0 if unknown) 
    P_cn_kW: Power losses [kW] in the core at nominal voltage
    P_const_kW: Power losses [kW] that is constant (bearing + friction and windage)
    """
    S_n_mva: float
    V_nom_kV: float
    cos_phi: float
    I_f_base: float
    R_a_nom: float
    R_f_nom: float
    X_d_u: float
    X_q_u: float
    X_l: float
    V_nom: float
    I_fd_pu_nom: float
    P_sn_kW: float
    P_rn_kW: float
    P_exn_kW: float
    P_cn_kW: float
    P_const_kW: float

    def __post_init__(self): 
        self.I_fd_pu_nom = self.I_fd_pu_nom
        self.P_sn = self.P_sn_kW / self.Sn_mva / 1000 # Convert to pu 
        self.P_rn = self.P_rn_kW / self.Sn_mva / 1000
        self.P_exn = self.P_exn_kW / self.Sn_mva / 1000
        self.P_cn = self.P_cn_kW / self.Sn_mva / 1000
        self.P_const = self.P_const_kW  / self.Sn_mva / 1000
        
        self.R_st = self.P_sn # Assumes I_a = 1.0 pu 
        self.R_rt = self.P_rn/self.I_fd_pu_nom**2 
      

class GeneratorModel2_If_in: 
    """ Main class for the generator loss model. Requires the model data. For loss calculations, field current need to be supplied. """
    def __init__(self, model_data: GenDataClass2) -> None: 
        self.md = model_data
        self.ex_states = ()
        self.alg_vars = ("P_g_pu", "Q_g_pu", "V_g_pu")

        P_g, Q_g, V_g, I_f = symbols("P_g Q_g V_g I_f")
        X_sym = (P_g, Q_g, V_g, I_f) 

        I_a_2 = (P_g**2 + Q_g**2)/V_g
        P_loss_st = I_a_2 * self.md.R_st
        P_loss_rt = self.md.R_rt * I_f**2
        P_loss_cr = self.md.P_cn*(V_g/self.md.V_nom)**2
        P_loss = P_loss_st + P_loss_rt + P_loss_cr + self.md.P_const

        dP_loss_st = [diff(P_loss_st, X) for X in X_sym] 
        dP_loss_rt = [diff(P_loss_rt, X) for X in X_sym]
        dP_loss_cr = [diff(P_loss_cr, X) for X in X_sym] 
        dP_loss = [P1+P2+P3 for P1, P2, P3 in zip(dP_loss_st, dP_loss_rt, dP_loss_cr)] 

        self.y = lambda X: np.array([lambdify(X_sym, P_loss)(*X)])
        self.dy_dx = lambda X: np.array([[lambdify(X_sym, dP_loss_i)(*X) for dP_loss_i in dP_loss]])

        self.P_loss_st = lambdify(X_sym, P_loss_st)
        self.P_loss_rt = lambdify(X_sym, P_loss_rt)
        self.P_loss_cr = lambdify(X_sym, P_loss_cr)

    
class GeneratorModel2: 
    """ Main class for the generator loss model. Requires the model data. For loss calculations, field current need to be supplied. """
    def __init__(self, model_data: GenDataClass2) -> None: 
        self.md = model_data
        self.ex_states = ()
        self.alg_vars = ("P_g_pu", "Q_g_pu", "V_g_pu")

        P_g, Q_g, V_g, I_f = symbols("P_g Q_g V_g I_f")
        X_sym = (P_g, Q_g, V_g, I_f) 

        I_a_2 = (P_g**2 + Q_g**2)/V_g
        P_loss_st = I_a_2 * self.md.R_st
        P_loss_rt = self.md.R_rt * I_f**2
        P_loss_cr = self.md.P_cn*(V_g/self.md.V_nom)**2
        P_loss = P_loss_st + P_loss_rt + P_loss_cr + self.md.P_const

        dP_loss_st = [diff(P_loss_st, X) for X in X_sym] 
        dP_loss_rt = [diff(P_loss_rt, X) for X in X_sym]
        dP_loss_cr = [diff(P_loss_cr, X) for X in X_sym] 
        dP_loss = [P1+P2+P3 for P1, P2, P3 in zip(dP_loss_st, dP_loss_rt, dP_loss_cr)] 

        self.y = lambda X: np.array([lambdify(X_sym, P_loss)(*X)])
        self.dy_dx = lambda X: np.array([[lambdify(X_sym, dP_loss_i)(*X) for dP_loss_i in dP_loss]])

        self.P_loss_st = lambdify(X_sym, P_loss_st)
        self.P_loss_rt = lambdify(X_sym, P_loss_rt)
        self.P_loss_cr = lambdify(X_sym, P_loss_cr)
