import numpy as np

class SatModelDataClass2: 
    """All values are in pu. """
    def __init__(self, bv: float, k: float, Cm: float, m: int, Ra_70: float, Xd_sat: float, Xq_sat: float, Xp: float) -> None: 
        self.bv = bv 
        self.k = k 
        self.Cm = Cm 
        self.m = m 
        self.Ra_70 = Ra_70 
        self.Xd_sat = Xd_sat 
        self.Xq_sat = Xq_sat 
        self.Xp = Xp 


class SaturationModel2: 
    """This saturation model is based on four parameters obtained from the open-circuit characteristic. Saturation is calculated through the Potier's reactance. """
    def __init__(self, model_data: SatModelDataClass2): 
        self.md = model_data
        self.ep_threshold = 0.55
        
    def calc_ifd(self, Vt, ia, phi) -> float: 
        delta = np.arctan(ia*(self.md.Xq_sat*np.cos(phi)-(self.md.Ra_70*np.sin(phi)))/(Vt+(self.md.Ra_70*ia*np.cos(phi))+self.md.Xq_sat*ia*np.sin(phi))) #Power load angle 
        egu = Vt*np.cos(delta) + (self.md.Ra_70*ia*np.cos(delta+phi)) + self.md.Xd_sat*ia*np.sin(delta+phi) 
        th = np.arctan(ia*(self.md.Xp*np.cos(phi) - self.md.Ra_70*np.sin(phi)) / (Vt + (self.md.Ra_70*ia*np.cos(phi)) + self.md.Xp*ia*np.sin(phi)) ) #Potiervinkel 
        ep = Vt*np.cos(th) + (self.md.Ra_70*ia*np.cos(th+phi)) + self.md.Xp*ia*np.sin(th+phi) 
        ifu = egu/self.md.bv #Field current at the air-gap line  
        ep_bool = ep > self.ep_threshold #Lager en array av True eller False (1 eller 0) slik at ifs kan lett regnes ut av i neste linje i mange tilfeller samtidig
        ifs = ((ep + self.md.Cm*ep**self.md.m)*self.md.k - ep/self.md.bv) * ep_bool #Field inkludert saturation
        ifd = (ifu + ifs) 
        return ifd, delta, egu