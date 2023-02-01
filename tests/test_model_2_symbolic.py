import unittest 
import numpy as np
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from SynGenLib.example_gens.model2 import Gen103MVA_sat_model_sym
from data_for_testing import data_103_MVA

class TestSatModel_symbolic(unittest.TestCase): 
    # test_gen = Gen103MVA
    # test_gen = Gen103MVA_If_in
    sat_model = Gen103MVA_sat_model_sym
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

    X_d_early = [0.8680258535554456, 0.8755206545714703, 0.8826507306986764, 0.8894200075994952, 0.8923902449680828, 0.8936652863395861,
                 0.8946631703926017, 0.8953854682379114]

    def test_nom_operation(self): 
        I_f_pred = self.sat_model.get_I_fd(self.P_tests[0], self.Q_tests[0], 1.0)
        self.assertAlmostEqual(I_f_pred*525.15, self.I_f_tests[0], places=-1)

    def test_all_points(self): 
        diff_meas = 0 
        diff_calc = 0
        for P, Q, If, If_calc in zip(self.P_tests, self.Q_tests, self.I_f_tests, self.I_f_calc_tests): 
            I_f_pred = self.sat_model.get_I_fd(P, Q, 1.0)
            diff_meas += (I_f_pred*525.15 - If)**2
            diff_calc += (I_f_pred*525.15 - If_calc)**2

        print(f"Sum of squared meas error: {diff_meas:.3f}, Sum of squared calc error: {diff_calc:.3f}")
        print(f"Yannick Calc model error: {442.320}")

    def test_y_selection(self): 
        #Compare X_d against early version. 
        for i, (P, Q) in enumerate(zip(self.P_tests, self.Q_tests)): 
            I_f, delta, X_d = self.sat_model.get_y_selection()(np.array([P, Q, 1.0])) # don't have any comparison here. 
            self.assertEqual(X_d, self.X_d_early[i])
    
if __name__ == "__main__": 
    unittest.main() 