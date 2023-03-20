import unittest 
import numpy as np
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from SynGenLib.example_gens.model3 import Gen103MVA_sat_model
from data_for_testing import data_103_MVA

class TestSatModel3(unittest.TestCase): 
    sat_model = Gen103MVA_sat_model
    P_tests =       data_103_MVA["P_tests"]
    cos_phi_tests = data_103_MVA["cos_phi_tests"]
    Q_tests =       data_103_MVA["Q_tests"]
    I_a_tests =     data_103_MVA["I_a_tests"]
    I_f_tests =     data_103_MVA["I_f_tests"]
    I_f_calc_tests =data_103_MVA["I_f_calc_tests"]
    P_f_tests =     data_103_MVA["P_f_tests"]
    P_ex_tests =    data_103_MVA["P_ex_tests"]
    P_br_tests =    data_103_MVA["P_br_tests"]
    P_a_tests =     data_103_MVA["P_a_tests"]
    P_s_tests =     data_103_MVA["P_s_tests"]
    P_c_tests =     data_103_MVA["P_c_tests"]
    P_be_tests =    data_103_MVA["P_be_tests"]
    P_wf_tests =    data_103_MVA["P_wf_tests"]
    P_loss_tests =  data_103_MVA["P_loss_tests"]
    eff_tests =     data_103_MVA["eff_tests"]
    I_f_base =      data_103_MVA["I_f_base"]
    I_a_base =      data_103_MVA["I_a_base"]
    phi_tests =     data_103_MVA["phi_tests"]

    def _test_sat_model(self): 
        Ifd_meas_error = []
        Ifd_calc_error = []
        for i in range(len(self.P_tests)): 
            I_fd, delta, psi_m = self.sat_model.calc_ifd(1.0, self.I_a_tests[i]/self.I_a_base, self.phi_tests[i])
            # print(f"Pred I_f: {I_fd*self.I_f_base}, meas = {self.I_f_tests[i]}")
            Ifd_meas_error.append((I_fd*self.I_f_base-self.I_f_tests[i])**2)
            Ifd_calc_error.append((I_fd*self.I_f_base-self.I_f_calc_tests[i])**2)

        print(f"Sum of squared meas error: {sum(Ifd_meas_error):.3f}, Sum of squared calc error: {sum(Ifd_calc_error):.3f}")
        print(f"Yannick Calc model error: {442.320}")

    def test_dany(self): 
        P = 93.562243/103
        Q = 47.669735/103
        V = 1.0
        I_a = np.sqrt(P**2 + Q**2)/V 
        phi = np.arctan(Q/P)
        ifd, delta, _ = self.sat_model.calc_ifd(V, I_a, phi)
        print(f"ifd={ifd*515.15}, delta = {delta}")

    
if __name__ == "__main__": 
    unittest.main() 