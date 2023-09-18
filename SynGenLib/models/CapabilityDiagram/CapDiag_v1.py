import numpy as np 
from numpy import cos, sin, sqrt, arctan
from scipy.interpolate import interp1d 
from scipy.optimize import root
from copy import deepcopy
from numba import njit 
import cmath as cm
from typing import Sequence, Tuple, Optional

from ..SatModels.SatModel2 import SaturationModel2

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.geometry.linestring import LineString

class CDFuncs1: 
    @staticmethod
    @njit
    def calc_Q1(V_g: float, xd_pu: float, xq_pu: float, E_min: float, P_g: float) -> float: 
        """Minimum internal voltage limit / UEL"""     
        c = (0.5*V_g**2 * (xd_pu+xq_pu))/(xd_pu*xq_pu) - E_min
        R = (0.5*V_g**2 * (xd_pu-xq_pu))/(xd_pu*xq_pu) + E_min 
        th = np.arcsin(P_g / R)
        Q_g = R*np.cos(th) - c
        return Q_g
    
    @staticmethod
    def calc_Q2(V_g: float, xd_pu: float, xq_pu: float, P_g: float, delta_max: float) -> float: 
        """Rotor Stability Limit"""
        def objective(X, delta_max): 
            Q_g, E_q = X
            a = E_q*xq_pu / (V_g*(xd_pu-xq_pu))
            b = np.sqrt(a**2 + 8) 
            cos_delta = np.clip(0.25 * (b - a), -1, 1) #stability factor ?
            delta = np.arccos(cos_delta)
            c = E_q*V_g / xd_pu * np.sin(delta) 
            d = 0.5*V_g**2 * (xd_pu - xq_pu)/xd_pu * xq_pu * np.sin(2*delta)
            P = c + d 
            Q = V_g*E_q/xd_pu*np.cos(delta_max) - V_g**2*(np.sin(delta_max)**2/xq_pu + np.cos(delta_max)**2/xd_pu)
            I_a = np.sqrt(P_g**2 + Q_g**2)/V_g
            phi = np.arctan(Q_g / P_g)    
            E = V_g + 1j*cm.rect(xd_pu*I_a, -phi) # + cm.rect(ra_pu*I_a, -phi)
            E = abs(E)
            f2 = Q - Q_g 
            f3 = E - E_q
            return np.array([f2, f3])

        sol = root(objective, x0=np.array([0.5, V_g]), args=(delta_max))
        Q_g, E_q = sol.x
        stab_m = 0.75 #Stability marigin (comes from??)
        return Q_g*stab_m
    
    @staticmethod
    @njit
    def calc_Q3(V_g: float, P_g: float, P_min: float, P_max: float) -> Tuple[float, float, float]:
        """Stator current limit \n 
        returns (P_new, Q_min, Q_max) \n 
        equations stems from I^2 = (P_g^2 + Q_g^2)/V_g^2 -> I = 1.0 to reach stator limit. \n 
        Therefore, P_g^2 + Q_g^2 = V_g^2 """
        if P_g < P_min:
            P_g = P_min 
        elif P_g > P_max: 
            P_g = P_max
        # P_g = np.clip(P_g, P_min, P_max)
        if P_g >= V_g: # Note that P_g cannot be bigger than V, even if Q = 0
            return (V_g, 0.0, 0.0)
        Q_max = np.sqrt(V_g**2 - P_g**2)
        Q_min = -Q_max # Because of symmetri
        return (P_g, Q_min, Q_max)
    
    @staticmethod   
    def calc_Q4(V_g: float, P_g: float, I_f_max: float, sat_model: SaturationModel2) -> float: 
        def objective(X): 
            Q = X[0] 
            I_a = sqrt(P_g**2 + Q**2)/V_g
            phi = arctan(Q/P_g)
            I_fd, delta, psi_m = sat_model.calc_ifd(V_g, I_a, phi)
            f1 = I_f_max - I_fd
            return np.array([f1])
        
        X0 = np.array([0.4])
        sol = root(objective, X0)
        Q_g = sol.x[-1]
        return Q_g

class CapabilityDataClass1: 
    """The data class given to the capability diagram model. """
    def __init__(self, P_lims: Tuple[float, float], X_d: float, X_q: float, R_a: float, I_f_max: Optional[float] = 2.0,  
                 E_min: Optional[float] = 0.1, delta_max: Optional[float] = 50*np.pi/180): 
        """P_lims: (P_min, P_max) in [pu] \n 
        X_d: (un)saturated synchronous reactance in d-axis [pu] \n 
        X_q: (un)saturated synchronous reactance in q-axis [pu] \n 
        R_a: Armature resistance at nominal operation [pu] \n 
        I_f_max: Maximum field current [pu] \n 
        E_min: Minimum internal voltage [pu] TODO: Change this? \n 
        delta_max: Maximum rotor angle before rotor stability threshold is reached [rad] """
        self.P_min = P_lims[0]
        self.P_max = P_lims[1]
        self.E_min = E_min if E_min is not None else 0.1
        self.I_f_max = I_f_max 
        self.delta_max = delta_max 
        self.X_d = X_d 
        self.X_q = X_q 
        self.R_a = R_a
    
    
class CapabilityDiagram1: 
    def __init__(self, cap_data: CapabilityDataClass1, sat_model: SaturationModel2): 
        self.md = cap_data
        self.sm = sat_model
    
    def get_Q_lims(self, V_g: float, P_pu: float) -> Tuple[float, float, float]: 
        """Finds the Q_pu limits (min, max) for a given voltage and active power \n 
        If the active power P_pu is outside the (P_min, P_max) range, a new P_pu will be returned which is \n
        inside the capability diagram. Respective Q_limits will also be returned.  
        Returns (P_new, Q_min, Q_max) """
        Q1 = CDFuncs1.calc_Q1(V_g, self.md.X_d, self.md.X_q, self.md.E_min, P_pu) 
        Q2 = CDFuncs1.calc_Q2(V_g, self.md.X_d, self.md.X_q, P_pu, self.md.delta_max) 
        P_g, Q3_min, Q3_max = CDFuncs1.calc_Q3(V_g, P_pu, self.md.P_min, self.md.P_max) 
        Q4 = CDFuncs1.calc_Q4(V_g, P_pu, self.md.I_f_max, self.sm)   
        Q_max = np.min((Q3_max, Q4))
        if np.isnan(Q1): 
            Q_min = np.max((Q2, Q3_min))
        else: 
            Q_min = np.max((Q1, Q2, Q3_min))
        return P_g, Q_min, Q_max 
    
    def is_inside(self, V: float, P_pu: float, Q_pu: float) -> bool: 
        P_new, Q_min, Q_max = self.get_Q_lims(V, P_pu)
        if P_pu > P_new or Q_min > Q_pu or Q_max < Q_pu: 
            return False 
        else: 
            return True


class CapabilityDiagram1Approx(CapabilityDiagram1): 
    """This class assumes no voltage dependency on the generator limits. Therefore, a polygon of 2*N_points is constructed using the module "shapely". \n 
    This is much faster compared to calculating the limits every time. """
    def __init__(self, cap_diag: CapabilityDiagram1, N_points: Optional[int]=100, V_nom: Optional[float]=1.0, tol: Optional[float]=1e-3): 
        """N_points: how many points that should make up the capability diagram. \n 
        V_nom: The nominal voltage used for making the limits of the capability diagram. \n
        tol: The allowable violation of Q limits."""
        self.CD = cap_diag 
        self.N_points = N_points
        self.V_nom = V_nom
        self.tol = tol
        points_1 = [] 
        points_2 = []
        for P in np.linspace(self.CD.md.P_min-2*self.tol, self.CD.md.P_max, N_points, endpoint=True):    
            P_new, Q_min, Q_max = self.CD.get_Q_lims(self.V_nom, P)
            # points_1.append(Point(Q_min, P_new))
            # points_2.append(Point(Q_max, P_new))
            points_1.append((Q_min - self.tol, P_new + self.tol))
            points_2.append((Q_max + self.tol, P_new + self.tol))
        self.CD_poly = Polygon(points_1 + points_2[::-1]) 
    
    def get_Q_lims(self, P_g: float) -> Tuple[float, float, float]:
        """Returns (P_new, Q_min, Q_max) """
        P_g = np.clip(P_g, self.CD.md.P_min, self.CD.md.P_max)
        line = LineString([Point(-2, P_g), Point(2, P_g)])
        lims = np.array(self.CD_poly.intersection(line).coords)
        if lims.shape[0] == 2:
            min_val, max_val = lims
            return P_g, min_val[0], max_val[0]
        else: 
            return P_g, 0.0, 0.0
    
    def is_inside(self, P_g, Q_g): 
        return self.CD_poly.contains(Point(Q_g, P_g))


