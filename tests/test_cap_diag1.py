import unittest 
import numpy as np
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from SynGenLib.example_gens.model2 import Gen103MVA_CD , Gen103MVA_CD_approx


class TestCD1(unittest.TestCase): 
    P_tests = [0.900, 0.675, 0.450, 0.225, 1.000, 0.750, 0.500, 0.250]
    cos_phi_tests = [0.9, 0.9, 0.9, 0.9, 1.0, 1.0, 1.0, 1.0]
    Q_tests = [np.tan(np.arccos(cos_phi))*P for cos_phi, P in zip(cos_phi_tests, P_tests)]  
    I_a_tests = [5406.1, 4054.6, 2703.0, 1351.5, 5406.1, 4054.6, 2703.0, 1351.5]
    I_f_tests = [1065.0, 936.12, 816.18, 711.38, 873.17, 776.61, 698.7, 646.84]
    I_f_calc_tests = [1064.88, 934.23, 815.01, 711.11, 857.50, 764.76, 691.88, 644.67]

    P_f_tests = [173.65, 133.66, 101.61, 77.19, 116.29, 91.99, 74.48, 63.81]
    P_ex_tests = [15.88, 13.02, 10.72, 8.87, 11.65, 9.92, 8.68, 7.92]
    P_br_tests = [2.13, 1.87, 1.63, 1.42, 1.75, 1.55, 1.40, 1.29]
    P_a_tests = [187.46, 105.45, 46.86, 11.72, 187.46, 105.45, 46.86, 11.72]
    P_s_tests = [89.16, 50.15, 22.30, 5.57, 89.16, 50.15, 22.30, 5.57]
    P_c_tests = [211.92]*8
    P_be_tests = [240.90]*8
    P_wf_tests = [172.92]*8
    P_loss_tests = [1094.02, 921.89, 808.85, 730.51, 1032.05, 884.81, 779.0, 172.92] # NOTE: Last index is wrong
    eff_tests = [98.834, 98.680, 98.250, 96.944, 99.008, 98.868, 98.509, 97.294]
    I_f_base = 525.15 # A

    def test_one_point_inside(self):
        self.assertTrue(Gen103MVA_CD.is_inside(1.0, 0.5, 0.0))

    def test_all_points_inside(self): 
        for i in range(len(self.P_tests)): 
            P_new, Q_min, Q_max = Gen103MVA_CD.get_Q_lims(1.0, self.P_tests[i])
            msg = f"Failed at P = {self.P_tests[i]}, Q = {self.Q_tests[i]}. Q_min = {Q_min}, Q_max = {Q_max}"
            self.assertTrue(Gen103MVA_CD.is_inside(1.0, self.P_tests[i], self.Q_tests[i]), msg=msg)

    def test_one_point_inside_approx(self): 
        self.assertTrue(Gen103MVA_CD_approx.is_inside(0.5, 0.0))

    def test_all_points_inside_approx(self): 
        for i in range(len(self.P_tests)): 
            P_new, Q_min, Q_max = Gen103MVA_CD_approx.get_Q_lims(self.P_tests[i])
            msg = f"Failed at P = {self.P_tests[i]}, Q = {self.Q_tests[i]}. Q_min = {Q_min}, Q_max = {Q_max}"
            self.assertTrue(Gen103MVA_CD_approx.is_inside(self.P_tests[i], self.Q_tests[i]), msg=msg)
    
if __name__ == "__main__": 
    unittest.main() 