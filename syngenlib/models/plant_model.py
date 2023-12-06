import numpy as np 
from typing import Sequence
import logging as lgn

from .generator_loss_model import *
from .transformer_loss_model import * 
from .capability_diagram_model import * 
from ..common import GeneratorDataClass, TransformerDataClass, CapabilityResult


class PlantModelResult: 
    pass 


class PlantModel: 
    """ Assumes each generator is connected through a separate transformer to a common busbar. """
    def __init__(self, 
                 gen_data: Sequence[GeneratorDataClass], 
                 trafo_data: Sequence[TransformerDataClass]): 
        if len(gen_data) != len(trafo_data): 
            lgn.warning("Different number of transformers and generators!")
            return 
        self.gen_data = gen_data
        self.trafo_data = trafo_data
        self.gen_models: Sequence[GeneratorLossModel] = []
        self.trafo_models: Sequence[TransformerLossModel] = []
        self.cap_diags: Sequence[CapabilityDiagram] = []

        # Change each trafo base to be the same as each generator base. 
        for gen, trafo in zip(self.gen_data, self.trafo_data): 
            trafo = trafo.change_base(gen.S_n_mva, gen.V_nom_kV, inplace=False) 
            self.gen_models.append(GeneratorLossModel(gen))
            self.trafo_models.append(TransformerLossModel(trafo))
            self.cap_diags.append(CapabilityDiagram(gen, trafo)) 

    def _get_Q_lims(self, P_s_mw: Sequence[Sequence[float]], V_hv: Sequence[float])->Sequence[CapabilityResult]: 
        """Given active power at the HV side of the generator transformers, 
        this calculates the reactive power limits. 
        Note: 
        P_s: Sequence[Sequence[float]] -> [P_g1_vals, P_g2_vals, ..., ]
        V_hv: Sequence[float] -> [V_hv_1, V_hv_2, ..., ]
        Returns (Q_min_tot_mvar, Q_max_tot_mvar)"""
        cd_res = []
        for P, cd, gen in zip(P_s_mw, self.cap_diags, self.gen_data): 
            cd_res.append(cd.calc_Q_lims(P/gen.S_n_mva, V_hv))
        return cd_res 
    
    def get_plant_Q_lims(self, P_s_mw: Sequence[Sequence[float]], V_hv: Sequence[float]): 
        """
        P_s_mw: Sequence[Sequence[float]] -> [P_g1_vals, P_g2_vals, ..., ]
        V_hv: Sequence[float] -> [V_hv_1, V_hv_2, ..., ]
        Returns (Q_min_tot_mvar, Q_max_tot_mvar)
        """
        cd_res = self._get_Q_lims(P_s_mw, V_hv) 
        Q_min_tot = 0 
        Q_max_tot = 0 
        for cd, gen in zip(cd_res, self.gen_data): 
            Q_min_tot += cd.Q_min_tot*gen.S_n_mva
            Q_max_tot += cd.Q_max_tot*gen.S_n_mva
        return Q_min_tot, Q_max_tot 
    
    def get_plant_losses(self, P_s_mw: Sequence[Sequence[float]], Q_s_mvar: Sequence[Sequence[float]], V_hv: Sequence[float]): 
        P_losses = []
        for P_s, Q_s, gen, trafo in zip(P_s_mw, Q_s_mvar, self.gen_models, self.trafo_models): 
            P = P_s/gen.md.S_n_mva 
            Q = Q_s/gen.md.S_n_mva
            I = (P - 1j*Q)/V_hv 
            V_g, I_g = trafo._calc_gen_V_I(V_hv, I) 
            S_g = V_g * I_g.conjugate()
            P_g = S_g.real
            Q_g = S_g.imag
            P_loss_trafo = trafo.get_P_losses_mw(P, Q, V_hv)
            P_loss_gen_res = gen.get_P_losses(P_g, Q_g, abs(V_g))
            print(P_loss_gen_res)
            P_loss_gen = P_loss_gen_res.get_losses_mw(gen.md.S_n_mva)
            P_loss_one = P_loss_trafo + P_loss_gen
            P_losses.append(P_loss_one)
        return P_losses 
