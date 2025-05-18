from ..data.data_types import GeneratorOperatingPoint, BranchOperatingPoint, PlantOperatingPoint
from ..data.result_classes import GeneratorBranchResults, PowerLossResult, CapabilityResults
from ..data.components import GeneratorDataclass, TransformerDataclass, GeneratorLossDataclass, CapabilityModelDataclass
from .saturation_model import SaturationBaseClass, LinearSaturationModel
from typing import Optional, Union
from math import sqrt, atan, atan2, nan, inf, cos, sin
import numpy as np 
from scipy.optimize import root 
from cmath import rect, phase

class GeneratorCalculationModel: 
    def __init__(self, gen_data: GeneratorDataclass, 
                 trafo_data: Optional[TransformerDataclass]=None, 
                 power_loss_data: Optional[GeneratorLossDataclass]=None,
                 saturation_model: Optional[SaturationBaseClass] = None, 
                 capability_model_data: Optional[CapabilityModelDataclass] = None):
        self.gen_data = gen_data
        P_nom_pu = gen_data.cos_phi 
        Q_nom_pu = sqrt(1.0 - P_nom_pu**2)
        nom_op = GeneratorOperatingPoint(P_mw=P_nom_pu*gen_data.S_n_mva,
                                         Q_mvar=Q_nom_pu*gen_data.S_n_mva,
                                         V_kv=gen_data.V_nom_kV)
        self.E_q_nom = self._get_E_q_from_gen_op(nom_op)

        if trafo_data is None: 
            self.trafo_data = TransformerDataclass.default_transformer(gen_data.S_n_mva, gen_data.V_nom_kV)
        else: 
            self.trafo_data = trafo_data
        
        if power_loss_data is None or not isinstance(power_loss_data, GeneratorLossDataclass): 
            self.power_loss_data = GeneratorLossDataclass.no_loss() 
        else: 
            self.power_loss_data = power_loss_data

        if saturation_model is None:            
            self.saturation_model = LinearSaturationModel(self.gen_data)
        else: 
            self.saturation_model = saturation_model 
        
        if capability_model_data is None:
            self.capability_data = CapabilityModelDataclass.default_limits(gen_data) 
        else: 
            self.capability_data = capability_model_data

        self.x_tot_pu = self.trafo_data.B.imag + self.gen_data.X_d_u 
        self._m = atan(self.capability_data.rotor_angle_max_rad)
        self._calculate_transformer_matrices()

    def change_tap_ratio(self, tap_ratio: float) -> None:
        self.trafo_data.change_tap_ratio(tap_ratio)
        self._calculate_transformer_matrices()

    def _calculate_transformer_matrices(self) -> None:
        S_base_change = (self.gen_data.S_n_mva/self.trafo_data.S_n_mva)
        V_base_change = (self.trafo_data.V_nom_lv_kV/self.gen_data.V_nom_kV)
        base_change = S_base_change * V_base_change**2
        self.Z_hv_pu = self.trafo_data.Z_hv * base_change
        self.Z_lv_pu = self.trafo_data.Z_lv * base_change
        self.Y_m_pu = self.trafo_data.Y_m / base_change

        self.A = 1.0 + self.Y_m_pu * self.Z_lv_pu 
        self.B = self.Z_lv_pu + self.Z_hv_pu + self.Y_m_pu * self.Z_lv_pu * self.Z_hv_pu 
        self.C = self.Y_m_pu 
        self.D = 1.0 + self.Y_m_pu * self.Z_hv_pu
        self.n_t = self.trafo_data.tap_ratio * self.gen_data.V_nom_kV / self.trafo_data.V_nom_lv_kV

        self.Y_bus = self.trafo_data.Y_bus / base_change 
        self.G_bus = self.Y_bus.real 
        self.B_bus = self.Y_bus.imag

    def get_E_q(self, op: Union[GeneratorOperatingPoint, BranchOperatingPoint, PlantOperatingPoint]) -> float:
        """Calculates the internal voltage E_q based on the generator's operating point.
        
        Args:
            op (Union[GeneratorOperatingPoint, BranchOperatingPoint, PlantOperatingPoint]):
                The operating point of the generator. 
            
        Returns:
            float: The internal voltage E_q [pu].
        """
        if type(op) == GeneratorOperatingPoint: 
            E_q = self._get_E_q_from_gen_op(op)
        elif type(op) == BranchOperatingPoint:
            E_q = self._get_E_q_from_branch_op(op)
        elif type(op) == PlantOperatingPoint:
            E_q = self._get_E_q_from_plant_op(op)
        else: 
            raise TypeError("operating_point must be of type GeneratorOperatingPoint, BranchOperatingPoint, or PlantOperatingPoint")
        return E_q

    def _get_E_q_from_gen_op(self, op: GeneratorOperatingPoint) -> float:
        P_g_pu = op.P_mw / self.gen_data.S_n_mva
        Q_g_pu = op.Q_mvar / self.gen_data.S_n_mva
        V_g_pu = op.V_kv / self.gen_data.V_nom_kV
        I_a = (P_g_pu - 1.0j*Q_g_pu) / V_g_pu
        E_q = abs(V_g_pu  + (self.gen_data.R_a + 1j*self.gen_data.X_q_u)*I_a)
        return E_q
    
    def _get_E_q_from_branch_op(self, op: BranchOperatingPoint) -> float:
        P_n_pu = op.P_mw / self.gen_data.S_n_mva
        Q_n_pu = op.Q_mvar / self.gen_data.S_n_mva
        V_n_pu = op.V_kv / self.trafo_data.V_nom_hv_kV 
        I_n = (P_n_pu - 1.0j*Q_n_pu) / V_n_pu 

        V_lv = V_n_pu / self.n_t 
        I_lv = I_n * self.n_t 
        V_g = V_lv * self.A + I_lv * self.B  
        I_g = V_lv * self.C + I_lv * self.D

        E_q = V_g + I_g*(self.gen_data.R_a + 1j*self.gen_data.X_q_u)
        E_q = abs(E_q)
        return E_q   
    
    def _get_E_q_from_plant_op(self, op: PlantOperatingPoint) -> float:
        branch_res = self._calculate_branch_results_from_plant_op(op) 
        gen_op = GeneratorOperatingPoint(branch_res.P_g_mw, branch_res.Q_g_mvar, branch_res.V_g_kv)
        E_q = self._get_E_q_from_gen_op(gen_op)
        return E_q
    
    def get_branch_results(self, operating_point: Union[GeneratorOperatingPoint, BranchOperatingPoint, PlantOperatingPoint]) -> GeneratorBranchResults: 
        if type(operating_point) == GeneratorOperatingPoint: 
            res = self._calculate_branch_results_from_gen_op(operating_point)
        elif type(operating_point) == BranchOperatingPoint:
            res = self._calculate_branch_results_from_branch_op(operating_point) 
        elif type(operating_point) == PlantOperatingPoint:
            res = self._calculate_branch_results_from_plant_op(operating_point)
        else: 
            raise TypeError("operating_point must be of type GeneratorOperatingPoint or BranchOperatingPoint")
        return res
        
    def _calculate_branch_results_from_gen_op(self, op: GeneratorOperatingPoint) -> GeneratorBranchResults: 
        n = self.trafo_data.tap_ratio

        V_n_gen = self.gen_data.V_nom_kV
        V_n_LV = self.trafo_data.V_nom_lv_kV 
        S_gen = self.gen_data.S_n_mva
        S_trafo = self.trafo_data.S_n_mva

        P_g = op.P_mw/self.gen_data.S_n_mva
        Q_g = op.Q_mvar/self.gen_data.S_n_mva
        V_g = op.V_kv/self.gen_data.V_nom_kV

        I_r = (P_g - 1.0j*Q_g)/V_g 
        I_lv = I_r * S_gen/S_trafo * V_n_LV/V_n_gen 
        V_lv = V_g * V_n_gen/V_n_LV 

        [V_n_t, I_n_t] = np.linalg.solve(self.trafo_data.M, [V_lv, I_lv])
        V_n = V_n_t * n 
        I_n = I_n_t / n 
        S_n = V_n * I_n.conjugate()
        P_n = S_n.real * self.trafo_data.S_n_mva
        Q_n = S_n.imag * self.trafo_data.S_n_mva

        return GeneratorBranchResults(op.P_mw, op.Q_mvar, P_n, Q_n, op.V_kv, abs(V_n)*self.trafo_data.V_nom_hv_kV)
    
    def _calculate_branch_results_from_branch_op(self, op: BranchOperatingPoint) -> GeneratorBranchResults: 
        n = self.trafo_data.tap_ratio
        V_n_gen = self.gen_data.V_nom_kV
        V_n_LV = self.trafo_data.V_nom_lv_kV 
        S_gen = self.gen_data.S_n_mva
        S_trafo = self.trafo_data.S_n_mva

        p_n = op.P_mw/self.trafo_data.S_n_mva 
        q_n = op.Q_mvar/self.trafo_data.S_n_mva
        V_n = op.V_kv/self.trafo_data.V_nom_hv_kV 

        i_hv = (p_n - 1.0j*q_n)/V_n 
        i_hv_t = i_hv * n 
        V_n_t = V_n / n 

        [V_lv, i_lv] = self.trafo_data.M @ np.array([V_n_t, i_hv_t])
        V_g = V_lv * V_n_LV/V_n_gen 
        i_g = i_lv * (S_trafo/S_gen) * (V_n_gen/V_n_LV) 
        s_g = V_g * i_g.conjugate()
        P_g = s_g.real * self.gen_data.S_n_mva
        Q_g = s_g.imag * self.gen_data.S_n_mva

        return GeneratorBranchResults(P_g, Q_g, op.P_mw, op.Q_mvar, abs(V_g)*self.gen_data.V_nom_kV, op.V_kv)
    
    def _calculate_branch_results_from_plant_op(self, op: PlantOperatingPoint) -> GeneratorBranchResults:
        V_n = op.V_n_kv / self.trafo_data.V_nom_hv_kV 
        V_lv = V_n / self.n_t 
        V_g_pu = op.V_g_kv / self.gen_data.V_nom_kV
        P_g_pu = op.P_mw / self.gen_data.S_n_mva

        def objective(delta_g): 
            P_g_calc_1 = V_g_pu**2*self.G_bus[0,0] 
            P_g_calc_2 = V_g_pu*V_lv*(self.G_bus[0,1]*cos(delta_g[0]) + self.B_bus[1,0]*sin(delta_g[0]))
            return (P_g_calc_1 + P_g_calc_2) - P_g_pu
        
        res = root(objective, np.array([0.0])) 
        if not res.success:
            raise ValueError("Root finding failed. Check the input values.")
        delta_g = res.x[0]
        V_g = rect(V_g_pu, delta_g) 
        v_vec = np.array([V_g, V_lv]) 
        i_inj = self.Y_bus @ v_vec 
        s_inj = v_vec * i_inj.conjugate() 
        P_g = s_inj[0].real * self.gen_data.S_n_mva
        Q_g = s_inj[0].imag * self.gen_data.S_n_mva
        P_n = -s_inj[1].real * self.gen_data.S_n_mva
        Q_n = -s_inj[1].imag * self.gen_data.S_n_mva
        return GeneratorBranchResults(P_g, Q_g, P_n, Q_n, op.V_g_kv, op.V_n_kv) 

    def get_branch_losses(self, op: Union[GeneratorOperatingPoint, BranchOperatingPoint, PlantOperatingPoint]) -> PowerLossResult: 
        if type(op) == GeneratorOperatingPoint: 
            res = self._calculate_branch_results_from_gen_op(op)
        elif type(op) == BranchOperatingPoint:
            res = self._calculate_branch_results_from_branch_op(op)
        elif type(op) == PlantOperatingPoint:
            res = self._calculate_branch_results_from_plant_op(op) 
        else: 
            raise TypeError("operating_point must be of type GeneratorOperatingPoint or BranchOperatingPoint")
        
        P_g = res.P_g_mw 
        Q_g = res.Q_g_mvar
        V_g = res.V_g_kv 
        gen_op = GeneratorOperatingPoint(P_g, Q_g, V_g)
        E_q_pu = self._get_E_q_from_gen_op(gen_op)
        I_f_pu = self.saturation_model.get_field_current(res, E_q_pu)
        I_g = sqrt(P_g**2 + Q_g**2) / V_g
        
        P_stator_loss_mw = self.power_loss_data.P_loss_nom_stator_pu * I_g**2 * self.gen_data.S_n_mva
        P_rotor_loss_mw = self.power_loss_data.P_loss_nom_rotor_pu * I_f_pu**2 * self.gen_data.S_n_mva
        P_core_loss_mw = self.power_loss_data.P_loss_nom_core_pu * V_g**2 * self.gen_data.S_n_mva
        P_const_loss_mw = self.power_loss_data.P_loss_nom_const_pu * self.gen_data.S_n_mva
        trafo_loss_mw = abs(res.P_g_mw - res.P_branch_mw)
        return PowerLossResult(P_stator_loss_mw, P_rotor_loss_mw, P_core_loss_mw, P_const_loss_mw, trafo_loss_mw)
    
    def get_Q_capability(self, operating_point: Union[GeneratorOperatingPoint, 
                                                      BranchOperatingPoint, 
                                                      PlantOperatingPoint]) -> CapabilityResults: 
        """Calculates the reactive power limits based on the generator's operational constraints.
        
        Args:
            operating_point (Union[GeneratorOperatingPoint, BranchOperatingPoint, PlantOperatingPoint]):
                The operating point of the generator. 
            
        Returns:
            CapabilityResults: A dataclass containing the reactive power limits and validation checks.

        Note: 
            You can not use the PlantOperatingPoint for this function.
        """
        if type(operating_point) == GeneratorOperatingPoint: 
            res = self._calculate_branch_results_from_gen_op(operating_point)
            P = res.P_branch_pu 
            V = res.V_grid_pu 

        elif type(operating_point) == BranchOperatingPoint: 
            P = operating_point.P_mw / self.gen_data.S_n_mva
            V = operating_point.V_pu 

        elif type(operating_point) == PlantOperatingPoint:
            raise ValueError("PlantOperatingPoint is not supported for this function.")
        
        else:
            raise TypeError("operating_point must be of type GeneratorOperatingPoint or BranchOperatingPoint")

        valid_P = self.capability_data.P_min_pu <= P <= self.capability_data.P_max_pu 
        Q_min_1, Q_max_1, valid_stator = self._get_stator_limits_pu(P, V) 
        Q_min_2, Q_max_2, valid_rotor = self._get_rotor_limits_pu(P, V)
        Q_min_3 = self._get_stability_limit_pu(P, V)
        Q_max_3 = inf # just for consistency with the comparisons. 
        Q_min_4, Q_max_4, valid_voltage = self._get_voltage_limits_pu(P, V)

        Q_min_vals = [Q_min_1, Q_min_2, Q_min_3, Q_min_4]
        Q_max_vals = [Q_max_1, Q_max_2, Q_max_3, Q_max_4]
        Q_min = max(Q_min_vals)
        Q_max = min(Q_max_vals)

        # The following code classifies which limiter is the most restrictive for the reactive power limits.
        limit_min = Q_min_vals.index(Q_min)
        limit_max = Q_max_vals.index(Q_max)

        cd_res = CapabilityResults(Q_min, Q_max, Q_min_1, Q_max_1, Q_min_2, Q_max_2, Q_min_3, Q_min_4, Q_max_4,
                                  valid_stator, valid_rotor, valid_P, valid_voltage, limit_min, limit_max)
        return cd_res

    def _get_rotor_limit_from_branch_op(self, op: BranchOperatingPoint) -> float:
        def obj(Q_n): 
            op_ = BranchOperatingPoint(op.P_mw, Q_n, op.V_kv) 
            E_q = self.get_E_q(op_)
            return E_q - self.capability_data.E_q_max 
        
        Q0 = sqrt(self.gen_data.S_n_mva** - op.P_mw**2)
        res = root(obj, np.array([Q0]))
        if not res.success:
            raise ValueError("Root finding for rotor limit failed.")
        Q_n = res.x[0]
        return Q_n
    
    def _get_rotor_limit_from_gen_op(self, op: GeneratorOperatingPoint) -> float:
        """TODO: This function 'cheats' by assuming no losses in the transformer. """
        branch_res = self.get_branch_results(op) 
        V_n = branch_res.V_grid_kv
        P_n = op.P_mw 
        op_ = BranchOperatingPoint(P_n, 0.0, V_n)
        Q_n = self._get_rotor_limit_from_branch_op(op_)
        return Q_n 
    
    def _get_rotor_limit_from_plant_op(self, op: PlantOperatingPoint) -> float:
        """TODO: This function 'cheats' by assuming no losses in the transformer. """
        V_n = op.V_n_kv
        P_n = op.P_mw 
        op_ = BranchOperatingPoint(P_n, 0.0, V_n)
        Q_n = self._get_rotor_limit_from_branch_op(op_)
        return Q_n
    
    def _get_rotor_limit(self, op: Union[GeneratorOperatingPoint, BranchOperatingPoint, PlantOperatingPoint]) -> float:
        """Calculates the maximum reactive power limit based on rotor current.

        Args:
            op (Union[GeneratorOperatingPoint, BranchOperatingPoint, PlantOperatingPoint]):
                The operating point.

        Returns:
            float: The maximum reactive power limited by the internal voltage [pu].
        """
        match op:
            case GeneratorOperatingPoint(): 
                Q_lim = self._get_rotor_limit_from_gen_op(op)
            case PlantOperatingPoint():
                Q_lim = self._get_rotor_limit_from_plant_op(op)
            case BranchOperatingPoint():
                Q_lim = self._get_rotor_limit_from_branch_op(op)
            case _:
                raise TypeError("operating_point must be of type GeneratorOperatingPoint, BranchOperatingPoint, or PlantOperatingPoint")
        return Q_lim

    def _get_stator_limits(self, op: Union[GeneratorOperatingPoint, BranchOperatingPoint, PlantOperatingPoint]) -> tuple[float, float, bool]:
        """Calculates the minimum and maximum reactive power limits based on stator current.
        
        Returns:
            tuple: A tuple containing (Q_min, Q_max, valid), where:
                Q_min (float): Minimum reactive power limit. [Mvar]
                Q_max (float): Maximum reactive power limit. [Mvar]
                valid (bool): Flag indicating the stator current limits are feasible.
        """

        match op: 
            case GeneratorOperatingPoint(): 
                Q_min, Q_max, valid_stator = self._get_stator_limits_from_gen_op(op)
            case PlantOperatingPoint():
                Q_min, Q_max, valid_stator = self._get_stator_limits_from_plant_op(op)
            case BranchOperatingPoint():
                Q_min, Q_max, valid_stator = self._get_stator_limits_from_branch_op(op)
            case _:
                raise TypeError("operating_point must be of type GeneratorOperatingPoint, BranchOperatingPoint, or PlantOperatingPoint") 
        if not valid_stator:
            raise ValueError("Stator limits are not valid. Active power is probably too high.")
        return Q_min, Q_max, valid_stator
        
    def _get_stator_limits_from_branch_op(self, op: BranchOperatingPoint) -> tuple[float, float, bool]:
        V_n_pu = op.V_kv / self.trafo_data.V_nom_hv_kV
        P_n_pu = op.P_mw / self.gen_data.S_n_mva 
        I_g_m = self.capability_data.I_a_max_pu
        
        a = (I_g_m/self.n_t * V_n_pu)**2 
        b = P_n_pu**2 
        if a < b: 
            valid_stator = False 
            return (nan, nan, valid_stator)
        Q_n_lim = sqrt(a - b)
        Q_min = -Q_n_lim*self.gen_data.S_n_mva
        Q_max = Q_n_lim*self.gen_data.S_n_mva
        valid_stator = True
        return (Q_min, Q_max, valid_stator)
    
    def _get_stator_limits_from_gen_op(self, op: GeneratorOperatingPoint) -> tuple[float, float, bool]:
        br_res = self.get_branch_results(op) 
        op_ = BranchOperatingPoint(br_res.P_branch_mw, 0.0, br_res.V_grid_kv)
        Q_min, Q_max, valid_stator = self._get_stator_limits_from_branch_op(op_)
        return (Q_min, Q_max, valid_stator)
    
    def _get_stator_limits_from_plant_op(self, op: PlantOperatingPoint) -> tuple[float, float, bool]:
        br_res = self.get_branch_results(op) 
        op_ = BranchOperatingPoint(br_res.P_branch_mw, 0.0, op.V_n_kv)
        Q_min, Q_max, valid_stator = self._get_stator_limits_from_branch_op(op_)
        return (Q_min, Q_max, valid_stator)
    
    def _get_stability_limit(self, op: Union[GeneratorOperatingPoint, BranchOperatingPoint, PlantOperatingPoint]) -> float:
        """Calculates the minimum reactive power limit based on maximum rotor angle.
        
        Returns:
            float: Minimum reactive power limit in Mvar [pu].
        """
        match op:
            case GeneratorOperatingPoint(): 
                Q_lim = self._get_stability_limit_from_gen_op(op)
            case PlantOperatingPoint():
                Q_lim = self._get_stability_limit_from_plant_op(op)
            case BranchOperatingPoint():
                Q_lim = self._get_stability_limit_from_branch_op(op)
            case _:
                raise TypeError("operating_point must be of type GeneratorOperatingPoint, BranchOperatingPoint, or PlantOperatingPoint")
        return Q_lim
    

    def _get_stability_limit_from_branch_op(self, op: BranchOperatingPoint) -> float:
        P_n_pu = op.P_mw / self.gen_data.S_n_mva
        V_n_pu = op.V_kv / self.trafo_data.V_nom_hv_kV
        V_lv = V_n_pu / self.n_t 
        X_tot = self.gen_data.X_d_u + self.B.imag
        Q_stab = self._m * P_n_pu - V_lv**2/X_tot 
        return Q_stab * self.gen_data.S_n_mva
    
    def _get_stability_limit_from_gen_op(self, op: GeneratorOperatingPoint) -> float:
        br_res = self.get_branch_results(op)
        op_ = BranchOperatingPoint(br_res.P_branch_mw, 0.0, br_res.V_grid_kv)
        Q_n = self._get_stability_limit_from_branch_op(op_)
        return Q_n
    
    def _get_stability_limit_from_plant_op(self, op: PlantOperatingPoint) -> float:
        br_res = self.get_branch_results(op)
        op_ = BranchOperatingPoint(br_res.P_branch_mw, 0.0, op.V_n_kv)
        Q_n = self._get_stability_limit_from_branch_op(op_)
        return Q_n
    
    def _get_voltage_limits(self, op) -> tuple[float, float, bool]:
        """Calculates the minimum and maximum reactive power limits based on voltage constraints.
        
        Returns:
            tuple: A tuple containing (Q_min_pu, Q_max_pu, valid), where:
                Q_min_pu (float): Minimum reactive power limit. [Mvar]
                Q_max_pu (float): Maximum reactive power limit. [Mvar]
                valid (bool): Flag indicating the voltage limits are feasible.
        """
        match op:
            case GeneratorOperatingPoint(): 
                Q_min, Q_max, valid = self._get_voltage_limits_from_gen_op(op)
            case PlantOperatingPoint():
                Q_min, Q_max, valid = self._get_voltage_limits_from_plant_op(op)
            case BranchOperatingPoint():
                Q_min, Q_max, valid = self._get_voltage_limits_from_branch_op(op)
            case _:
                raise TypeError("operating_point must be of type GeneratorOperatingPoint, BranchOperatingPoint, or PlantOperatingPoint")
        return (Q_min, Q_max, valid)
            
    def _get_voltage_limits_from_branch_op(self, op: BranchOperatingPoint) -> tuple[float, float, bool]:
        V_n_pu = op.V_kv / self.trafo_data.V_nom_hv_kV
        V_lv = V_n_pu / self.n_t
        P_n_pu = op.P_mw / self.gen_data.S_n_mva 
        V_min = self.capability_data.V_g_min
        V_max = self.capability_data.V_g_max
        X_T = self.B.imag 
        
        d_min = np.arcsin(P_n_pu * X_T / (V_min*V_lv))
        d_max = np.arcsin(P_n_pu * X_T / (V_max*V_lv))
        Q_n_min = V_min*V_lv*cos(d_min)/(X_T) - V_lv**2/(X_T)
        Q_n_max = V_max*V_lv*cos(d_max)/(X_T) - V_lv**2/(X_T)
        return (Q_n_min*self.gen_data.S_n_mva, Q_n_max*self.gen_data.S_n_mva, True)


    def _get_voltage_limits_from_gen_op(self, op: GeneratorOperatingPoint) -> tuple[float, float, bool]:
        br_res = self.get_branch_results(op)
        op_ = BranchOperatingPoint(br_res.P_branch_mw, 0.0, br_res.V_grid_kv)
        Q_min, Q_max, valid = self._get_voltage_limits_from_branch_op(op_)
        return (Q_min, Q_max, valid)

    def _get_voltage_limits_from_plant_op(self, op: PlantOperatingPoint) -> tuple[float, float, bool]:
        br_res = self.get_branch_results(op)
        op_ = BranchOperatingPoint(br_res.P_branch_mw, 0.0, op.V_n_kv)
        Q_min, Q_max, valid = self._get_voltage_limits_from_branch_op(op_)
        return (Q_min, Q_max, valid)
