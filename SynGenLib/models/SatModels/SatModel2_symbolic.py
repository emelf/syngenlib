import numpy as np
from numpy import cos, sin
from typing import Tuple, Sequence
from scipy.optimize import root
from dataclasses import dataclass 
import sympy as sp 
from sympy import diff, lambdify, symbols 

from .SatModel2 import SatModelDataClass2


class SaturationModel2_sym: 
    """
    Class for definint the saturation equations of the synchronous generator. 
    
    Attributes 
    ----------
    model_data: SatModelDataClass2

    Naming convention 
    ----------
    u: Inputs = [P_g, Q_g, V_g]
    y: All algebraic variables used for calculations, given in the y_vars attribute
    y_selection: Algebraic variables relevant for this model. Chosen to be [I_f, delta, X_d]. 

    Methods
    ----------
    get_y: Returns a lambda function where the input is a Sequence of inputs, and returns a numpy array of all the states. TODO - Implement this. 
    get_y_selection: Returns a lambda function of y_selection with one input -> u: Sequence[P_g, Q_g, V_g] 
    get_dy_du_selection: Returns a numpy matrix/array of the partial derivatives of y w.r.t. u -> shape = (N_y, N_u) 
    get_I_fd: If only field current calculation is required. 
    get_dIf_dQ: If only partial derivative of If w.r.t. Q_g is required. 
    """
    def __init__(self, model_data: SatModelDataClass2): 
        self.md = model_data
        self.x_vars = {}
        self.y_vars = {"I_a": "Armature current", "phi": "Phase angle between voltage and current", 
                       "delta": "Rotor angle", "e_d": "d-axis terminal voltage", "e_q": "q-axis terminal voltage", 
                       "i_d": "d-axis armature current", "i_q": "q-axis armature current", "psi_d": "d-axis flux linkage", 
                       "psi_q": "q-axis flux linkage", "psi_m": "flux linkage magnitude", "c_sat": "Saturation coefficient", 
                       "sat_d": "saturation coefficient in the d-axis", "X_ad": "Saturated mutual reactance between d-axis and armature", 
                       "X_d": "saturated d-axis reactance", "I_f": "Field current in the rotor circuit"}
        self.y_selection = ("I_f", "delta", "X_d")
        self.u_vars = {"P_g": "Generator active power [pu]", "Q_g": "Generator reactive power [pu]", "V_g": "Generator terminal voltage"}

        self.dy_du_vars =(["dIf_dPg", "dIf_dQg", "dIf_dVg"], 
                          ["ddelta_dPg", "ddelta_dQg", "ddelta_dVg"], 
                          ["dX_d_sat_dPg", "dX_d_sat_dQg", "dX_d_sat_dVg"])

        self._generate_y_equations() 

    def _generate_y_equations(self): 
        P_g, Q_g, V_g = sp.symbols("P_g Q_g V_g")
        self._u_symbols = (P_g, Q_g, V_g)

        self.phi = sp.atan(Q_g/P_g)
        self.I_a = sp.sqrt(P_g**2 + Q_g**2)/V_g
        num = self.I_a*(self.md.X_q_u*sp.cos(self.phi)-(self.md.R_a*sp.sin(self.phi)))
        denum = (V_g+(self.md.R_a*self.I_a*sp.cos(self.phi))+self.md.X_q_u*self.I_a*sp.sin(self.phi))
        self.delta = sp.atan(num/denum)

        self.e_d = V_g * sp.sin(self.delta)
        self.e_q = V_g * sp.cos(self.delta) 
        self.i_d = self.I_a * sp.sin(self.delta + self.phi) 
        self.i_q = self.I_a * sp.cos(self.delta + self.phi) 
        self.psi_d = self.e_q + self.md.R_a*self.i_q 
        self.psi_q = -self.e_d - self.md.R_a*self.i_d 
        self.psi_m = sp.sqrt((self.psi_d + self.md.X_l*self.i_d)**2 + (self.psi_q + self.md.X_l*self.i_q)**2)

        self.c_sat = self.md.SG10 * self.psi_m**self.md.exp / self.psi_m
        self.sat_d = 1 / (1 + self.c_sat)

        self.X_ad_sat = self.md.X_ad_u * self.sat_d
        self.X_d_sat = self.X_ad_sat + self.md.X_l 
        self.I_fd = (self.e_q + self.md.R_a*self.i_q + self.X_d_sat*self.i_d) / self.X_ad_sat 

        dIf_dQ = diff(self.I_fd, Q_g) 
        self.get_I_fd = lambdify(self._u_symbols, self.I_fd)
        self.get_dIf_dQ = lambdify(self._u_symbols, dIf_dQ)

    def get_y(self): 
        """
        Obtains all equations for y, with u as input. \n
        In this case, u = (P_g, Q_g, V_g)
        """
        #TODO: Finish this 
        pass 

    def get_y_selection(self) -> Sequence[float]: 
        y = lambda u: np.array([lambdify(self._u_symbols, self.I_fd)(*u), 
                                lambdify(self._u_symbols, self.delta)(*u), 
                                lambdify(self._u_symbols, self.X_d_sat)(*u)])
        return y 

    def get_dy_du_selection(self): 
        dIf_dX = [diff(self.I_fd, u) for u in self._u_symbols]
        ddelta_dX = [diff(self.delta, u) for u in self._u_symbols]
        dXd_dX = [diff(self.X_d_sat, u) for u in self._u_symbols]
        dy_du = lambda u: np.array([[lambdify(self._u_symbols, dX_dY)(*u) for dX_dY in dIf_dX], 
                                    [lambdify(self._u_symbols, dX_dY)(*u) for dX_dY in ddelta_dX],
                                    [lambdify(self._u_symbols, dX_dY)(*u) for dX_dY in dXd_dX]])
        return dy_du


