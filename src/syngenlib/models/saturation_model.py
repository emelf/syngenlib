from abc import ABC, abstractmethod
from math import sqrt, atan2, cos, sin
from ..data.components import GeneratorDataclass
from ..data.result_classes import GeneratorBranchResults 

class SaturationBaseClass(ABC):
    """
    Abstract base class for saturation models.
    """

    def __init__(self, *args, **kwargs):
        pass # Initialize any common attributes or methods here

    @abstractmethod
    def get_field_current(self, P_g_pu: float, Q_g_pu: float, V_g_pu: float, gen_data: GeneratorDataclass) -> float:
        """ 
        Returns the per-unit field current based on the generator's active power, reactive power, and voltage."""
        return NotImplementedError("This method should be overridden by subclasses")
    

class LinearSaturationModel(SaturationBaseClass): 
    """
    Linear saturation model for synchronous generators.

    This model assumes a linear relation between the internal voltage E_q and the field current I_f. 
    The calculation of the field current is I_f = k_f * E_q, where k_f is a constant. 

    Attributes:
        k_f (float): Field current constant [pu]
    """

    def __init__(self, gen_data: GeneratorDataclass): 
        super().__init__()
        P_g_pu = gen_data.cos_phi
        Q_g_pu = sqrt(1 - P_g_pu**2) 
        V_g_pu = 1.0 
        # E_q = self._get_E_q(P_g_pu, Q_g_pu, V_g_pu, gen_data) TODO 
        E_q = 1.6 
        self.k_If = 1.0/E_q # Assumes nominal operating point has I_f = 1.0 pu. 
        self.S_base_mva = gen_data.S_n_mva

    def get_field_current(self, gen_res: GeneratorBranchResults) -> float:
        I_f = self.k_If * gen_res.E_q_pu
        return I_f