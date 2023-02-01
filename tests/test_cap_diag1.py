import unittest 
import numpy as np
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from SynGenLib.example_gens.model2 import Gen103MVA_CD , Gen103MVA_CD_approx
from data_for_testing import data_103_MVA

class TestCD1(unittest.TestCase): 
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