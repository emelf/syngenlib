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
            trafo = trafo.change_base(gen.S_n_mva, inplace=False) 
            self.gen_models.append(GeneratorLossModel(gen))
            self.trafo_models.append(TransformerDataClass(trafo))
            self.cap_diags.append(CapabilityDiagram(gen, trafo)) 

    def get_Q_lims(self, P_s: Sequence[float], V_hv): 
        """Given active power at the HV side of the generator transformers, 
        this calculates the reactive power limits. 
        Returns (Q_min_tot_mvar, Q_max_tot_mvar)"""
        cd_res = []
        for P, cd in zip(P_s, self.cap_diags): 
            cd_res.append(cd.calc_Q_lims(P, V_hv))
        return cd_res 

    

    
      