import numpy as np
from numpy import cos, sin
from typing import Tuple
from scipy.optimize import root

class SatModelDataClass2: 
    """All values are in pu. """
    def __init__(self, X_d_u: float, X_q_u: float, X_l: float, R_a: float, SG10: float, SG12: float) -> None: 
        """X_d_u -> Unsaturated d-axis reactance [pu] \n 
        X_l -> Leakage reactance [pu] \n 
        R_a -> Nominal stator resistance [pu] \n 
        SG10 -> Saturation coefficient at nominal voltage [.] \n 
        SG12 -> Saturation coefficient at 1.2 pu voltage [.] \n """
        self.X_d_u = X_d_u
        self.X_q_u = X_q_u
        self.X_l = X_l
        self.R_a = R_a 
        self.SG10 = SG10 
        self.SG12 = SG12 

        self.X_ad_u = self.X_d_u - self.X_l
        self.exp = np.log(1.2 * self.SG12 / self.SG10) / np.log(1.2)


# class SaturationModel2: 
#     """NOTE: Assumes the power factory model for saturation (exponential model) and that the q-axis reactance remains unsaturated.  """
#     def __init__(self, model_data: SatModelDataClass2): 
#         self.md = model_data
        
#     def calc_ifd(self, V_t, I_a, phi) -> Tuple[float, float, float]: 
#         """returns (I_fd, delta, psi_m)"""
#         delta = np.arctan(I_a*(self.md.X_q_u*np.cos(phi)-(self.md.R_a*np.sin(phi)))/(V_t+(self.md.R_a*I_a*np.cos(phi))+self.md.X_q_u*I_a*np.sin(phi))) #Power load angle 
#         e_d = V_t * sin(delta)
#         e_q = V_t * cos(delta) 
#         i_d = I_a * sin(delta + phi) 
#         i_q = I_a * cos(delta + phi) 
#         psi_d = e_q + self.md.R_a*i_q 
#         psi_q = -e_d - self.md.R_a*i_d 
#         psi_m = np.sqrt((psi_d + self.md.X_l*i_d)**2 + (psi_q + self.md.X_l*i_q)**2)

#         c_sat = self.md.SG10 * psi_m**self.md.exp / psi_m
#         sat_d = 1 / (1 + c_sat)
#         X_ad_sat = self.md.X_ad_u * sat_d
#         X_d_sat = X_ad_sat + self.md.X_l 
#         # print(f"X_d_sat = {X_d_sat}")
#         I_fd = (e_q + self.md.R_a*i_q + X_d_sat*i_d) / X_ad_sat 
#         return I_fd, delta, psi_m

class SaturationModel2: 
    """NOTE: Assumes the power factory model for saturation (exponential model) and that the q-axis reactance remains unsaturated.  """
    def __init__(self, model_data: SatModelDataClass2): 
        self.md = model_data

    def _objective(self, X, V_t, I_a, phi):
        delta, c_sat, sat_d = X 
        e_d = V_t * sin(delta)
        e_q = V_t * cos(delta) 
        i_d = I_a * sin(delta + phi) 
        i_q = I_a * cos(delta + phi) 
        psi_d = e_q + self.md.R_a*i_q 
        psi_q = -e_d - self.md.R_a*i_d 
        psi_m = np.sqrt((psi_d + self.md.X_l*i_d)**2 + (psi_q + self.md.X_l*i_q)**2)

        f1 = delta - np.arctan(I_a*(self.md.X_q_u*np.cos(phi)-(self.md.R_a*np.sin(phi)))/(V_t+(self.md.R_a*I_a*np.cos(phi))+self.md.X_q_u*I_a*np.sin(phi)))
        f2 = c_sat - self.md.SG10 * psi_m**self.md.exp / psi_m
        f3 = sat_d - 1 / (1 + c_sat)
        return np.array([f1, f2, f3])
        
    def calc_ifd(self, V_t, I_a, phi) -> Tuple[float, float, float]: 
        """returns (I_fd, delta, psi_m)"""
        #Initial calculation
        delta = np.arctan(I_a*(self.md.X_q_u*np.cos(phi)-(self.md.R_a*np.sin(phi)))/(V_t+(self.md.R_a*I_a*np.cos(phi))+self.md.X_q_u*I_a*np.sin(phi)))
        e_d = V_t * sin(delta)
        e_q = V_t * cos(delta) 
        i_d = I_a * sin(delta + phi) 
        i_q = I_a * cos(delta + phi) 
        psi_d = e_q + self.md.R_a*i_q 
        psi_q = -e_d - self.md.R_a*i_d 
        psi_m = np.sqrt((psi_d + self.md.X_l*i_d)**2 + (psi_q + self.md.X_l*i_q)**2)
        
        c_sat = self.md.SG10 * psi_m**self.md.exp / psi_m
        sat_d = 1 / (1 + c_sat)

        delta_res, c_sat_res, sat_d_res = root(self._objective, np.array([delta, c_sat, sat_d]), args=(V_t, I_a, phi)).x

        X_ad_sat = self.md.X_ad_u * sat_d_res
        X_d_sat = X_ad_sat + self.md.X_l 
        # print(f"X_d_sat = {X_d_sat}")
        I_fd = (e_q + self.md.R_a*i_q + X_d_sat*i_d) / X_ad_sat 
        return I_fd, delta_res, psi_m