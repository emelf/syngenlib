from typing import Optional
import numpy as np

from ..common import TransformerDataClass, TransformerLossResult

class TransformerLossModel: 
    def __init__(self, trafo_data: TransformerDataClass): 
        self.md = trafo_data 

    def _calc_gen_V_I(self, V_r, I_r): 
        VI = np.array([V_r, I_r])
        [V_g, I_g] = self.md.ABCD_mat @ VI 
        return V_g, I_g 

    def _calc_P_losses(self, V_r, I_r, V_s, I_s): 
        I_12 = I_r + V_r*self.md.Y_hv 
        P_hv_loss = abs(V_r)**2 * self.md.Y_hv.real 
        P_lv_loss = abs(V_s)**2 * self.md.Y_lv.real 
        P_12_loss = abs(I_12)**2 * self.md.Z_12.real 
        return P_hv_loss + P_lv_loss + P_12_loss 
    
    def get_P_losses(self, P_r_pu, Q_r_pu, V_r) -> TransformerLossResult: 
        """All input quantities are refered to at the grid-side (recieving side)"""
        I_grid = (P_r_pu - 1j*Q_r_pu)/V_r
        V_g, I_g = self._calc_gen_V_I(V_r, I_grid) 
        # P_loss = self._calc_P_losses(V_r, I_grid, V_g, I_g)
        S_gen = V_g * I_g.conj() 
        return TransformerLossResult(S_gen.real, P_r_pu, S_gen.imag, Q_r_pu)

        

        