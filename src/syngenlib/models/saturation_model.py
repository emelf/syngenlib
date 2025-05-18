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
    def get_field_current(self, branch_res: GeneratorBranchResults, E_q: float) -> float:
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
        Q_g_pu = sqrt(1.0 - P_g_pu**2)
        V_g_pu = 1.0 
        I_a = (P_g_pu - 1.0j*Q_g_pu) / V_g_pu
        self.E_q_nom = abs(V_g_pu  + (gen_data.R_a + 1j*gen_data.X_q_u)*I_a)
        self.k_If = 1.0/self.E_q_nom  

    def get_field_current(self, branch_res: GeneratorBranchResults, E_q_pu: float) -> float:
        I_f_pu = self.k_If * E_q_pu
        return I_f_pu