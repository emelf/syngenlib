import unittest 
import numpy as np
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from SynGenLib.example_gens.model2 import Gen103MVA_sat_model_sym


class TestSatModel_symbolic(unittest.TestCase): 
    # test_gen = Gen103MVA
    # test_gen = Gen103MVA_If_in
    sat_model = Gen103MVA_sat_model_sym
    P_tests = [0.900, 0.675, 0.450, 0.225, 1.000, 0.750, 0.500, 0.250]
    cos_phi_tests = [0.9, 0.9, 0.9, 0.9, 1.0, 1.0, 1.0, 1.0]
    Q_tests = [np.tan(np.arccos(cos_phi))*P for cos_phi, P in zip(cos_phi_tests, P_tests)]  
    phi_tests = [np.arctan(Q/P) for Q, P in zip(Q_tests, P_tests)]
    I_a_tests = [5406.1, 4054.6, 2703.0, 1351.5, 5406.1, 4054.6, 2703.0, 1351.5]
    I_f_tests = [1065.0, 936.12, 816.18, 711.38, 873.17, 776.61, 698.7, 646.84]
    I_f_tests = [1055.0, 936.12, 816.18, 711.38, 873.17, 776.61, 698.7, 646.84]
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
    I_f_base = 525.15
    I_a_base = 5406

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