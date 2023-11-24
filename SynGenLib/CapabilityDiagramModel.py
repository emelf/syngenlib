import numpy as np 
from typing import Tuple, Optional
from .DataClasses import GenDataClass, TrafoDataClass

class CapabilityDiagram: 
    """This is the generator+trafo capability diagram constructor. 
    It may also work as the generator capability diagram only when no trafo 
    model is supplied. When defining the terminal voltage, this could be either 
    the generator voltage, or the upstream trafo voltage, depending on if the 
    trafo model has been given or not. """
    def __init__(self, gen_data: GenDataClass, trafo_data: Optional[TrafoDataClass]=None): 
        self.gen_data = gen_data    
        if trafo_data is None: 
            self.trafo_data = TrafoDataClass(S_n_mva=self.gen_data.S_n_mva, 
                                             V_nom_kV=self.gen_data.V_nom_kV, 
                                             V_SCH=0.0, I_E=0.0, P_Cu=0.0, P_Fe=0.0)
        else: 
            self.trafo_data = trafo_data
        self.m = np.arctan(self.gen_data.delta_max) 
        self.X_tot = self.trafo_data.X_T + self.gen_data.X_d_u

    def _calc_stator_limit(self, P, V) -> Tuple[float, float]:
        Q_max = np.sqrt((V*self.gen_data.I_g_max/self.trafo_data.tap_ratio)**2 - P**2)
        return (-Q_max, Q_max)        
    
    def _calc_rotor_limit(self, P, V) -> Tuple[float, float]: 
        r_f = self.gen_data.E_q_max*V/(self.trafo_data.tap_ratio*self.X_tot)
        q_f = -V**2/(self.trafo_data.tap_ratio**2 * self.X_tot)
        Q_g_max = np.sqrt(r_f**2 - P**2) + q_f 
        Q_g_min = -np.sqrt(r_f**2 - P**2) + q_f 
        return (Q_g_min, Q_g_max) 
    
    def _calc_stab_limit(self, P, V) -> Tuple[float, float]: 
        c = -V**2/(self.trafo_data.tap_ratio**2*self.X_tot)
        Q_min = self.m * P + c
        return (Q_min, None)
    
    def _calc_voltage_limit(self, P, V):
        if self.trafo_data.X_T <= 0: 
            return (None, None)
        k1 = V/(self.trafo_data.tap_ratio*self.trafo_data.X_T)
        Q_max = k1*(self.gen_data.V_g_max - V/self.trafo_data.tap_ratio)
        Q_min = k1*(self.gen_data.V_g_min - V/self.trafo_data.tap_ratio)
        return (Q_min, Q_max)
    
    def calc_Q_lims(self, P_g, V_g) -> Tuple[float, float]: 
        """ Returns (Q_min, Q_max, valid). 
        The valid flag is a boolean signal indicating if the 
        simulation is valid or not. An invalid simulation may be 
        when the active power limit is exceeded. """
        is_valid = self.gen_data.P_g_min_pu <= P_g <= self.gen_data.P_g_max_pu
        Q_min_1, Q_max_1 = self._calc_stator_limit(P_g, V_g)
        _, Q_max_2 = self._calc_rotor_limit(P_g, V_g)
        Q_min_3, _ = self._calc_stab_limit(P_g, V_g)
        Q_min = max((Q_min_1, Q_min_3))
        Q_max = min((Q_max_1, Q_max_2))
        return (Q_min, Q_max, is_valid)
    

if __name__ == "__main__": 
    import matplotlib.pyplot as plt 
    from ExampleGens import gen_103_mva
    
    CD1 = CapabilityDiagram(gen_103_mva)

    P_vals = np.linspace(1e-3, 1, 100)
    V_g = 1.0 
    Q_min_stator = []
    Q_max_stator = []
    Q_min_rotor = []
    Q_max_rotor = []
    Q_min_stab = [] 
    Q_max_stab = []
    Q_min_v = []
    Q_max_v = [] 
    Q_min_vals = [] 
    Q_max_vals = []
    
    for P in P_vals: 
        Q_min_1, Q_max_1 = CD1._calc_stator_limit(P, V_g)
        Q_min_stator.append(Q_min_1)
        Q_max_stator.append(Q_max_1)

        Q_min_2, Q_max_2 = CD1._calc_rotor_limit(P, V_g)
        Q_max_rotor.append(Q_max_2)

        Q_min_3, Q_max_3 = CD1._calc_stab_limit(P, V_g)
        Q_min_stab.append(Q_min_3)

        Q_min_4, Q_max_4 = CD1._calc_voltage_limit(P, V_g) 
        Q_min_v.append(Q_min_4)
        Q_max_v.append(Q_max_4)

        Q_min_vals.append(max(Q_min_1, Q_min_3))
        Q_max_vals.append(min(Q_max_1, Q_max_2))

    fig1 = plt.figure(1)
    plt.plot(Q_min_stator, P_vals, color="black", label="Stator limit")
    plt.plot(Q_max_stator, P_vals, color="black")
    plt.plot(Q_max_rotor, P_vals, color="orange", label="Rotor limit")
    plt.plot(Q_min_stab, P_vals, color="blue", label="Angle limit")
    plt.plot(Q_min_v, P_vals, color="grey", linestyle="dashed", label="V_min")
    plt.plot(Q_max_v, P_vals, color="grey", linestyle="dashed", label="V_max")

    plt.legend()
    plt.grid()

    fig2 = plt.figure(2)
    plt.plot(Q_min_vals, P_vals, color="black")
    plt.plot(Q_max_vals, P_vals, color="black")
    plt.grid()

    plt.show()

    
        
