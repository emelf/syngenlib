import numpy as np 

class CapDiag: 
    def __init__(self, X_d, E_q_max, delta_max, P_g_min, P_g_max): 
        self.X_d = X_d 
        self.E_q_max = E_q_max
        self.delta_max = delta_max
        self.P_g_min = P_g_min 
        self.P_g_max = P_g_max

        self.r_f_1 = self.E_q_max/self.X_d 
        self.q_f_1 = -1.0/self.X_d
        self.m = np.arctan(self.delta_max) 

    def calc_stator_limit(self, P_g, V_g) -> (float, float):
        if P_g <= V_g: 
            Q_max = np.sqrt(1.0*V_g - P_g**2) 
            return (-Q_max, Q_max)
        elif P_g > self.P_g_max:
            return (0, 0)
        
    
    def calc_rotor_limit(self, P_g, V_g) -> (float, float): 
        r_f = self.r_f_1 * V_g 
        q_f = self.q_f_1*V_g**2 
        Q_g_max = np.sqrt(r_f**2 - P_g**2) + q_f 
        Q_g_min = -np.sqrt(r_f**2 - P_g**2) + q_f 
        return (Q_g_min, Q_g_max) 
    
    def calc_stab_limit(self, P_g, V_g) -> (float, float): 
        Q_g_min = self.m * P_g + self.q_f_1*V_g**2
        return (Q_g_min, None)
    

if __name__ == "__main__": 
    import matplotlib.pyplot as plt 
    CD1 = CapDiag(1.0, 1.8, np.pi/180*30, 0.1, 0.9)
    P_vals = np.linspace(1e-3, 1, 100)
    V_g = 1.0 
    Q_min_stator = []
    Q_max_stator = []
    Q_min_rotor = []
    Q_max_rotor = []
    Q_min_stab = [] 
    Q_max_stab = [] 
    Q_min_vals = [] 
    Q_max_vals = []
    
    for P in P_vals: 
        Q_min_1, Q_max_1 = CD1.calc_stator_limit(P, V_g)
        Q_min_stator.append(Q_min_1)
        Q_max_stator.append(Q_max_1)

        Q_min_2, Q_max_2 = CD1.calc_rotor_limit(P, V_g)
        Q_max_rotor.append(Q_max_2)

        Q_min_3, Q_max_3 = CD1.calc_stab_limit(P, V_g)
        Q_min_stab.append(Q_min_3)

        Q_min_vals.append(max(Q_min_1, Q_min_3))
        Q_max_vals.append(min(Q_max_1, Q_max_2))

    fig1 = plt.figure(1)
    plt.plot(Q_min_stator, P_vals, color="black", label="Stator limit")
    plt.plot(Q_max_stator, P_vals, color="black")
    plt.plot(Q_max_rotor, P_vals, color="orange", label="Rotor limit")
    plt.plot(Q_min_stab, P_vals, color="blue", label="Angle limit")
    plt.legend()
    plt.grid()

    fig2 = plt.figure(2)
    plt.plot(Q_min_vals, P_vals, color="black")
    plt.plot(Q_max_vals, P_vals, color="black")
    plt.grid()

    plt.show()

    
        
