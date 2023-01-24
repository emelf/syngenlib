import unittest 
import numpy as np
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from SynGenLoss.GenModel2.example_gens import Gen103MVA 


class TestModel2(unittest.TestCase): 
    test_gen = Gen103MVA
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

    def test_nom_loss(self):
        res_loss = self.test_gen.get_P_losses(self.P_tests[0], self.Q_tests[0], 1.0)

        self.assertAlmostEqual(res_loss.P_core_pu*self.test_gen.md.Sn_mva*1000, 
                               self.P_c_tests[0], places=1)

        self.assertAlmostEqual(res_loss.P_rotor_pu*self.test_gen.md.Sn_mva*1000, 
                               self.P_ex_tests[0] + self.P_br_tests[0] + self.P_f_tests[0], places=1)

        self.assertAlmostEqual(res_loss.P_stator_pu*self.test_gen.md.Sn_mva*1000, 
                               self.P_a_tests[0] + self.P_s_tests[0], places=1)

        self.assertAlmostEqual(res_loss.P_const_pu*self.test_gen.md.Sn_mva*1000, 
                               self.P_be_tests[0] + self.P_wf_tests[0], places=1) 

    def test_all_abjora(self): 
        for i in range(len(self.P_tests)): 
            res_loss = self.test_gen.get_P_losses(self.P_tests[i], self.Q_tests[i], 1.0)

            self.assertAlmostEqual(res_loss.P_core_pu*self.test_gen.md.Sn_mva*1000, 
                                   self.P_c_tests[i], 
                                   places=-1, msg=f"Error at idx={i}")

            self.assertAlmostEqual(res_loss.P_rotor_pu*self.test_gen.md.Sn_mva*1000, 
                                   self.P_ex_tests[i] + self.P_br_tests[i] + self.P_f_tests[i], 
                                   places=-1, msg=f"Error at idx={i}")

            self.assertAlmostEqual(res_loss.P_stator_pu*self.test_gen.md.Sn_mva*1000, 
                                   self.P_a_tests[i] + self.P_s_tests[i], 
                                   places=-1, msg=f"Error at idx={i}")

            self.assertAlmostEqual(res_loss.P_const_pu*self.test_gen.md.Sn_mva*1000, 
                                   self.P_be_tests[i] + self.P_wf_tests[i], 
                                   places=-1, msg=f"Error at idx={i}") 

    def test_calc_currents(self): 
        for i in range(len(self.P_tests)): 
            ia, ifd, delta = self.test_gen._calc_currents(self.P_tests[i], self.Q_tests[i], 1.0)
            
            self.assertAlmostEqual(ia * self.test_gen.md.Ia_nom_A, self.I_a_tests[i], places=1, msg=f"I_a error at idx={i}")
            self.assertAlmostEqual(ifd * self.test_gen.md.I_f_base, self.I_f_calc_tests[i], places=0, msg=f"I_f error at idx={i}")
    
    
if __name__ == "__main__": 
    unittest.main() 