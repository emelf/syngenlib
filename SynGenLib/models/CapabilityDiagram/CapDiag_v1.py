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
            cos_delta = 0.25 * (b - a) #stability factor ?
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
    def __init__(self, cap_diag: CapabilityDiagram1, N_points: Optional[int]=100, V_nom: Optional[float]=1.0): 
        self.CD = cap_diag 
        self.N_points = N_points
        self.V_nom = V_nom
        points_1 = [] 
        points_2 = []
        for P in np.linspace(self.CD.md.P_min, self.CD.md.P_max, N_points):    
            P_new, Q_min, Q_max = self.CD.get_Q_lims(self.V_nom, P)
            points_1.append(Point(Q_min, P_new))
            points_2.append(Point(Q_max, P_new))
        self.CD_poly = Polygon(points_1 + points_2[::-1]) 
    
    def get_Q_lims(self, P_g: float) -> Tuple[float, float, float]:
        """Returns (P_new, Q_min, Q_max) """
        P_g = np.clip(P_g, self.CD.md.P_min, self.CD.md.P_max)
        line = LineString([Point(-2, P_g), Point(2, P_g)])
        lims = np.array(self.CD_poly.intersection(line).coords)
        if lims.shape[0] == 2:
            min_val, max_val = lims
            return min_val[0], max_val[0]
        else: 
            return 0.0, 0.0
    
    def is_inside(self, P_g, Q_g): 
        return self.CD_poly.contains(Point(Q_g, P_g))


# class TrafoCapabilityDiagram(CapabilityDiagram): 
#     def __init__(self, gen: GeneratorModel, trafo: TrafoModel): 
#         super().__init__(gen)
#         self.trafo = deepcopy(trafo)
#         self.trafo_md = deepcopy(trafo.md)
#         self.X_T_pu = self.trafo_md.X_T
        
#     def get_Q_lims_plant(self, V_hv: float, P_pu: float): 
#         """Return -> (P_g, Q_hv_min, Q_hv_max, V_g_min, V_g_max)"""
#         def objective(X, min_max, P_g_in): 
#             V_g, Q_hv_lim = X 
#             P_g, Q_g_min, Q_g_max = self.get_Q_lims(V_g, P_g_in)
#             if min_max == 0: 
#                 Q_g_lim = Q_g_min 
#             else: 
#                 Q_g_lim = Q_g_max 
#             delta_g = self._calc_delta_g(V_hv, V_g, P_g)
#             Q_g = self._calc_Q_g(V_hv, V_g, delta_g)
#             Q_hv = self._calc_Q_hv(V_g, P_g, Q_g_lim)
#             f1 = Q_g_lim - Q_g
#             f2 = Q_hv_lim - Q_hv 
#             return np.array([f1, f2])
        
#         P_g, Q_min_0, Q_max_0 = self.get_Q_lims(V_hv, P_pu)
#         X0_min = np.array([V_hv, Q_min_0])
#         X0_max = np.array([V_hv, Q_max_0]) 
#         sol_min = root(objective, X0_min, args=(0, P_pu)) 
#         sol_max = root(objective, X0_max, args=(1, P_pu)) 
#         V_g_min, Q_hv_min = sol_min.x 
#         V_g_max, Q_hv_max = sol_max.x 
        
#         P_g_min, Q_g_min, _ = self.get_Q_lims(V_g_min, P_pu)
#         P_g_max, _, Q_g_max = self.get_Q_lims(V_g_max, P_pu)
#         return P_g_min, P_g_max, Q_g_min, Q_g_max, Q_hv_min, Q_hv_max, V_g_min, V_g_max
        
#     def calc_Q_at_Vg_lim(self, V_hv, V_g_lim, P_g): 
#         delta = self._calc_delta_g(V_hv, V_g_lim, P_g)
#         Q_g = self._calc_Q_g(V_hv, V_g_lim, delta)
#         Q_g_hv = self._calc_Q_hv(V_g_lim, P_g, Q_g)
#         return Q_g_hv
    
#     def is_inside(self, V: float, P_pu: float, Q_pu: float, V_min=0.95, V_max=1.05) -> bool: 
#         P = P_pu * self.md.Sn_mva 
#         Q = Q_pu * self.md.Sn_mva
#         P_g, Q_g, V_g, _ = self.trafo.calc_PQV_sending(P, Q, V)
#         P_g = P_g / self.md.Sn_mva
#         Q_g = Q_g / self.md.Sn_mva
#         P_new, Q_min, Q_max = self.get_Q_lims(V_g, P_g)
#         P_g_cond = P_g <= P_new
#         Q_g_cond = Q_min <= Q_g <= Q_max
#         V_g_cond =  V_min <= V_g <= V_max
#         return P_g_cond and Q_g_cond and V_g_cond

#     def _calc_Q_g(self, V_hv, V_g, delta): 
#         return V_g**2/self.X_T_pu - V_hv*V_g*cos(delta)/self.X_T_pu 
    
#     def _calc_delta_g(self, V_hv, V_g, P_g): 
#         return np.arcsin(P_g * self.X_T_pu / (V_g * V_hv)) 
    
#     def _calc_Q_hv(self, V_g, P_g, Q_g): 
#         return (Q_g*V_g - self.X_T_pu*(P_g**2 + Q_g**2)) / (V_g**2) 
    
  
# class PlantCapabilityDiagram: 
#     def __init__(self, gens: Sequence[GeneratorModel], trafos: Sequence[TrafoModel]): 
#         self.gen_mds = [deepcopy(gen.md) for gen in gens] 
#         self.sat_mds = [deepcopy(gen.satmodel) for gen in gens]
#         self.trafo_mds = [deepcopy(trafo.md) for trafo in trafos]
#         self.X_Ts_pu = np.array([trafo_md.X_T for trafo_md in self.trafo_mds])
#         self.CDs = [CapabilityDiagram(gen) for gen in gens] 
#         self.CDs_trafo = [TrafoCapabilityDiagram(gen, trafo) for gen, trafo in zip(gens, trafos)]
#         self.Sn_base = np.array([md.Sn_mva for md in self.gen_mds])
        
#     def get_Q_lims_plant(self, V_hv: float, P_gs_mva: Sequence[float]): 
#         """Return -> (P_g, Q_hv_min, Q_hv_max, V_g_min, V_g_max) \n 
#         Note: Assuming equal pu distribution. """
#         Q_gs_min =  np.zeros(len(self.gen_mds))
#         Q_gs_max =  np.zeros(len(self.gen_mds))
#         Q_hvs_min = np.zeros(len(self.gen_mds))
#         Q_hvs_max = np.zeros(len(self.gen_mds))
#         V_gs_min =  np.zeros(len(self.gen_mds))
#         V_gs_max =  np.zeros(len(self.gen_mds))
#         P_gs_min =  np.zeros(len(self.gen_mds))
#         P_gs_max =  np.zeros(len(self.gen_mds))
        
#         def assign_to(idx, P_g_min, P_g_max, Q_g_min, Q_g_max, Q_hv_min, Q_hv_max, V_g_min, V_g_max): 
#             Q_gs_min[idx] = Q_g_min
#             Q_gs_max[idx] = Q_g_max
#             Q_hvs_min[idx] = Q_hv_min
#             Q_hvs_max[idx] = Q_hv_max
#             V_gs_min[idx] = V_g_min
#             V_gs_max[idx] = V_g_max
#             P_gs_min[idx] = P_g_min
#             P_gs_max[idx] = P_g_max
        
#         for i, (CD_trafo, P_hv_mva) in enumerate(zip(self.CDs_trafo, P_gs_mva)): 
#             S_base = CD_trafo.md.Sn_mva
#             P_pu = P_hv_mva / S_base 
#             P_g_min, P_g_max, Q_g_min, Q_g_max, Q_hv_min, Q_hv_max, V_g_min, V_g_max = CD_trafo.get_Q_lims_plant(V_hv, P_pu)
#             assign_to(i, P_g_min*S_base, P_g_max*S_base, Q_g_min*S_base, Q_g_max*S_base, 
#                       Q_hv_min*S_base, Q_hv_max*S_base, V_g_min, V_g_max)
            
#         return P_gs_min, P_gs_max, Q_gs_min, Q_gs_max, Q_hvs_min, Q_hvs_max, V_gs_min, V_gs_max
    
#     def get_Q_plant(self, V_hv: float, V_gs: Sequence[float], P_gs: Sequence[float]): 
#         """Input: V_hv, V_gs -> Return (Q_plant, Q_gs, Q_gs_hv) \n 
#         V_hv -> The voltage at the HV busbar \n 
#         V_gs -> A sequence of the generator terminal voltages. """
#         P_gs = P_gs / self.Sn_base
#         deltas = np.array([self._calc_delta_g(V_hv, V_g, P_g, X_T) for V_g, P_g, X_T in zip(V_gs, P_gs, self.X_Ts_pu)] )
#         Q_gs = np.array([self._calc_Q_g(V_hv, V_g, delta, X_T) for V_g, delta, X_T in zip(V_gs, deltas, self.X_Ts_pu)] ) * self.Sn_base
#         Q_hvs = np.array([self._calc_Q_hv(V_g, P_g, Q_g, X_T) for V_g, P_g, Q_g, X_T in zip(V_gs, P_gs, Q_gs, self.X_Ts_pu)] ) * self.Sn_base
#         Q_plant = Q_hvs.sum() 
#         return Q_plant, Q_gs, Q_hvs
        
#     def calc_Q_at_Vg_lim(self, V_hv, V_g_lims: Sequence[float], P_gs: Sequence[float]) -> Sequence[float]: 
#         deltas = np.array([self._calc_delta_g(V_hv, V_g, P_g, X_T) for V_g, P_g, X_T in zip(V_g_lims, P_gs, self.X_Ts_pu)] )
#         Q_gs = np.array([self._calc_Q_g(V_hv, V_g, delta, X_T) for V_g, delta, X_T in zip(V_g_lims, deltas, self.X_Ts_pu)] ) * self.Sn_base
#         Q_hvs = np.array([self._calc_Q_hv(V_g, P_g, Q_g, X_T) for V_g, P_g, Q_g, X_T in zip(V_g_lims, P_gs, Q_gs, self.X_Ts_pu)] ) * self.Sn_base
#         return Q_hvs
    
#     def is_inside(self, V_hv: float, P_hvs_mva: Sequence[float], Q_hvs_mva: Sequence[float]) -> bool: 
#         P_hvs_pu = np.array(P_hvs_mva) / self.Sn_base
#         Q_hvs_pu = np.array(Q_hvs_mva) / self.Sn_base
#         G1_inside = self.CDs_trafo[0].is_inside(V_hv, P_hvs_pu[0], Q_hvs_pu[0])
#         G2_inside = self.CDs_trafo[1].is_inside(V_hv, P_hvs_pu[1], Q_hvs_pu[1])
#         return G1_inside & G2_inside
    
#     def _calc_delta_g(self, V_hv, V_g, P_g, X_T): 
#         return np.arcsin(P_g * X_T / (V_g * V_hv)) 
        
#     def _calc_Q_g(self, V_hv, V_g, delta, X_T): 
#         return V_g**2/X_T - V_hv*V_g*cos(delta)/X_T

#     def _calc_Q_hv(self, V_g, P_g, Q_g, X_T): 
#         return (Q_g*V_g - X_T*(P_g**2 + Q_g**2)) / (V_g**2) 
