class GenLossRes: 
    def __init__(self, P_in: float, P_rotor: float, P_ex: float, P_stator: float, P_core: float, P_const: float) -> None: 
        self.P_in = P_in
        self.P_const_pu = P_const 
        self.P_rotor_pu = P_rotor 
        self.P_ex_pu = P_ex
        self.P_stator_pu = P_stator 
        self.P_core_pu = P_core 
        self.P_loss_tot = P_rotor + P_ex + P_stator + P_core + P_const
        self.eff = P_in / (P_in + self.P_loss_tot) 
        
    def __str__(self): 
        str1 = f"P_in = {self.P_in:.3f}, P_loss_tot = {self.P_loss_tot:.3f}, eff = {self.eff*100:.3f} %"
        str2 = f"P_l_stator = {self.P_stator_pu:.3f}, P_l_rotor = {self.P_rotor_pu:.3f}, P_l_ex = {self.P_ex_pu:.3f}, P_l_core = {self.P_core_pu:.3f}, P_const = {self.P_const_pu:.3f}"
        return str1 + "\n" + str2

    def __mul__(self, S_base: float): 
        new_genlossres = GenLossRes(self.P_in*S_base, self.P_rotor_pu*S_base, self.P_ex_pu*S_base, self.P_stator_pu*S_base, 
                                    self.P_core_pu*S_base, self.P_const_pu*S_base)
        return new_genlossres
