import numpy as np 
from typing import Tuple, Optional
from .DataClasses import GenDataClass, TrafoDataClass, CapabilityResult, CapabilityLimit

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
            self.no_trafo = True # Flag that tells that no trafo should be considered 
        else: 
            self.trafo_data = trafo_data
            self.no_trafo = False 
        self.m = np.arctan(self.gen_data.delta_max) 
        self.X_tot = self.trafo_data.X_T + self.gen_data.X_d_u

    def _calc_stator_limit(self, P, V) -> Tuple[float, float, bool]:
        Q_max = np.zeros_like(P, dtype=float)
        valid_stator = P**2 <= (V*self.gen_data.I_g_max/self.trafo_data.tap_ratio)**2
        Q_max[valid_stator] = np.sqrt((V[valid_stator]*self.gen_data.I_g_max/self.trafo_data.tap_ratio)**2 - P[valid_stator]**2)
        return (-Q_max, Q_max, valid_stator)   
    
    def _calc_rotor_limit(self, P, V) -> Tuple[float, float]: 
        r_f_max = self.gen_data.E_q_max*V/(self.trafo_data.tap_ratio*self.X_tot)
        r_f_min = self.gen_data.E_q_min*V/(self.trafo_data.tap_ratio*self.X_tot)
        q_f = -V**2/(self.trafo_data.tap_ratio**2 * self.X_tot)
        below_min = P < r_f_min
        valid = P <= r_f_max
        valid_and_below_min = np.logical_and(below_min, valid)
        Q_g_max = np.zeros_like(P, dtype=float)
        Q_g_min = np.array([np.nan for _ in Q_g_max], dtype=float)

        Q_g_max[valid] = np.sqrt(r_f_max[valid]**2 - P[valid]**2) + q_f[valid] 
        Q_g_min[valid_and_below_min] = np.sqrt(r_f_min[valid_and_below_min]**2 - P[valid_and_below_min]**2) + q_f[valid_and_below_min]
        return (Q_g_min, Q_g_max, valid) 
    
    def _calc_stab_limit(self, P, V) -> Tuple[float, float]: 
        c = -V**2/(self.trafo_data.tap_ratio**2*self.X_tot)
        Q_min = self.m * P + c
        Q_max = np.array([np.nan for _ in Q_min], dtype=float)
        return (Q_min, Q_max)
    
    def _calc_voltage_limit(self, P, V):
        Q_min = np.zeros_like(P, dtype=float)
        Q_max = np.zeros_like(P, dtype=float)
        if self.trafo_data.X_T <= 0: 
            Q_min[:] = np.nan 
            Q_max[:] = np.nan
        k1 = V/self.trafo_data.tap_ratio/self.trafo_data.X_T 
        k2_min = (self.gen_data.V_g_min*k1)
        k2_max = (self.gen_data.V_g_max*k1)
        valid = k2_min >= P 
        Q_min[valid] = np.sqrt(k2_min**2 - P**2) - k1*V/self.trafo_data.tap_ratio
        Q_max[valid] = np.sqrt(k2_max**2 - P**2) - k1*V/self.trafo_data.tap_ratio
        Q_min[~valid] = np.nan
        Q_max[~valid] = np.nan
        return (Q_min, Q_max, valid)
    
    def calc_Q_lims(self, P_g, V_g) -> CapabilityResult: 
        """ Returns res: CapabilityResult """
        valid_power = np.logical_and(self.gen_data.P_g_min_pu <= P_g, 
                                     P_g <= self.gen_data.P_g_max_pu) 
        Q_min_1, Q_max_1, valid_stator = self._calc_stator_limit(P_g, V_g)
        Q_min_2, Q_max_2, valid_rotor = self._calc_rotor_limit(P_g, V_g)
        Q_min_3, Q_max_3 = self._calc_stab_limit(P_g, V_g)

        if self.no_trafo:
            Q_min_4 = np.array([np.nan for _ in range(len(valid_power))]) 
            Q_max_4 = np.array([np.nan for _ in range(len(valid_power))]) 
            valid_voltage = np.logical_and(self.gen_data.V_g_min <= V_g, 
                                           self.gen_data.V_g_max >= V_g)
        else: 
            Q_min_4, Q_max_4, valid_voltage = self._calc_voltage_limit(P_g, V_g)

        limit_min = np.nanargmax([Q_min_1, Q_min_2, Q_min_3, Q_min_4], axis=0)
        limit_max = np.nanargmin([Q_max_1, Q_max_2, Q_max_3, Q_max_4], axis=0)
        Q_min = np.nanmax([Q_min_1, Q_min_2, Q_min_3, Q_min_4], axis=0)
        Q_max = np.nanmin([Q_max_1, Q_max_2, Q_max_3, Q_max_4], axis=0)
        valid_voltage = np.logical_and.reduce((Q_min_4 <= Q_max, Q_max_4 >= Q_min, valid_voltage))

        res = CapabilityResult(P_g, Q_min, Q_max, Q_min_1, Q_max_1, Q_min_2, Q_max_2, Q_min_3, Q_min_4, Q_max_4, 
                               valid_stator, valid_rotor, valid_power, valid_voltage, 
                               limit_min, limit_max)
        return res
    
    def calc_Q_lims_mesh(self, p_mesh, v_mesh) -> CapabilityResult: 
        Q_min = np.zeros_like(p_mesh, dtype=float) 
        Q_max = np.zeros_like(p_mesh, dtype=float) 
        Q_min_1 = np.zeros_like(p_mesh, dtype=float) 
        Q_max_1 = np.zeros_like(p_mesh, dtype=float) 
        Q_min_2 = np.zeros_like(p_mesh, dtype=float) 
        Q_max_2 = np.zeros_like(p_mesh, dtype=float) 
        Q_min_3 = np.zeros_like(p_mesh, dtype=float) 
        Q_min_4 = np.zeros_like(p_mesh, dtype=float) 
        Q_max_4 = np.zeros_like(p_mesh, dtype=float) 
        valid_stator = np.zeros_like(p_mesh, dtype=bool) 
        valid_rotor = np.zeros_like(p_mesh, dtype=bool) 
        valid_power = np.zeros_like(p_mesh, dtype=bool) 
        valid_voltage = np.zeros_like(p_mesh, dtype=bool) 
        limit_min = np.zeros_like(p_mesh, dtype=int) 
        limit_max = np.zeros_like(p_mesh, dtype=int)
        for i, (p, v) in enumerate(zip(p_mesh, v_mesh)): 
            res = self.calc_Q_lims(p, v)
            Q_min[i] = res.Q_min_tot
            Q_max[i] = res.Q_max_tot
            Q_min_1[i] = res.Q_stator_min
            Q_max_1[i] = res.Q_stator_max
            Q_min_2[i] = res.Q_rotor_min
            Q_max_2[i] = res.Q_rotor_max
            Q_min_3[i] = res.Q_stab_min
            Q_min_4[i] = res.Q_v_min
            Q_max_4[i] = res.Q_v_max
            valid_stator[i] = res.valid_stator_current
            valid_rotor[i] = res.valid_rotor_current
            valid_power[i] = res.valid_active_power
            valid_voltage[i] = res.valid_voltage_levels
            limit_min[i] = res.limiter_Q_min
            limit_max[i] = res.limiter_Q_max
        return CapabilityResult(p_mesh, Q_min, Q_max, Q_min_1, Q_max_1, Q_min_2, Q_max_2, Q_min_3, Q_min_4, Q_max_4, 
                                valid_stator, valid_rotor, valid_power, valid_voltage, 
                                limit_min, limit_max) 
    

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

    
        
