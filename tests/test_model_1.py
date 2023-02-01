import unittest 
import numpy as np
import numdifftools as nd 
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from SynGenLib.example_gens.model1 import Gen103MVA 
from data_for_testing import data_103_MVA

class TestModel1(unittest.TestCase): 
    test_gen = Gen103MVA
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

    def test_nom_loss(self):
        res_loss = self.test_gen.get_P_losses(self.P_tests[0], self.Q_tests[0], 1.0)

        self.assertAlmostEqual(res_loss.P_core_pu*self.test_gen.md.Sn_mva*1000, 
                               self.P_c_tests[0], places=1)

        self.assertAlmostEqual((res_loss.P_rotor_pu+res_loss.P_ex_pu)*self.test_gen.md.Sn_mva*1000, 
                               self.P_ex_tests[0] + self.P_br_tests[0] + self.P_f_tests[0], places=1)

        self.assertAlmostEqual(res_loss.P_stator_pu*self.test_gen.md.Sn_mva*1000, 
                               self.P_a_tests[0] + self.P_s_tests[0], places=1)

        self.assertAlmostEqual(res_loss.P_const_pu*self.test_gen.md.Sn_mva*1000, 
                               self.P_be_tests[0] + self.P_wf_tests[0], places=1) 

    def test_all_abjora(self): 
        calc_eff = []
        for i in range(len(self.P_tests)): 
            res_loss = self.test_gen.get_P_losses(self.P_tests[i], self.Q_tests[i], 1.0)
            calc_eff.append(res_loss.eff*100)
            self.assertAlmostEqual(res_loss.P_core_pu*self.test_gen.md.Sn_mva*1000, 
                                   self.P_c_tests[i], 
                                   places=-1, msg=f"Error at idx={i}")

            self.assertAlmostEqual((res_loss.P_rotor_pu+res_loss.P_ex_pu)*self.test_gen.md.Sn_mva*1000, 
                                   self.P_ex_tests[i] + self.P_br_tests[i] + self.P_f_tests[i], 
                                   places=-1, msg=f"Error at idx={i}")

            self.assertAlmostEqual(res_loss.P_stator_pu*self.test_gen.md.Sn_mva*1000, 
                                   self.P_a_tests[i] + self.P_s_tests[i], 
                                   places=-1, msg=f"Error at idx={i}")

            self.assertAlmostEqual(res_loss.P_const_pu*self.test_gen.md.Sn_mva*1000, 
                                   self.P_be_tests[i] + self.P_wf_tests[i], 
                                   places=-1, msg=f"Error at idx={i}") 
        calc_eff = np.array(calc_eff) 
        eff_diff = calc_eff - np.array(self.eff_tests)
        print(f"Efficiency abs diff = {np.sum(np.abs(eff_diff))}")

    def test_calc_currents(self): 
        for i in range(len(self.P_tests)): 
            ia, ifd, delta = self.test_gen._calc_currents(self.P_tests[i], self.Q_tests[i], 1.0)
            
            self.assertAlmostEqual(ia * self.test_gen.md.Ia_nom_A, self.I_a_tests[i], places=1, msg=f"I_a error at idx={i}")
            self.assertAlmostEqual(ifd * self.test_gen.md.If_nom_A, self.I_f_calc_tests[i], places=0, msg=f"I_f error at idx={i}")

    def test_sat_model(self): 
        Ifd_meas_error = []
        Ifd_calc_error = []
        for i in range(len(self.P_tests)): 
            ia, ifd, delta = self.test_gen._calc_currents(self.P_tests[i], self.Q_tests[i], 1.0)
            Ifd_meas_error.append((ifd*self.I_f_base-self.I_f_tests[i])**2)
            Ifd_calc_error.append((ifd*self.I_f_base-self.I_f_calc_tests[i])**2)

        print(f"Sum of squared meas error: {sum(Ifd_meas_error):.3f}, Sum of squared calc error: {sum(Ifd_calc_error):.3f}")

    def test_dany(self): 
        ia, ifd, delta = Gen103MVA._calc_currents(0.8985, 0.2583, 1.0)
        print(f"ia = {ia}, ifd={ifd*525.15}, delta = {delta}")


def gen_loss_wrapper_Q_g(P_g: float, V_g: float, gen):
     return lambda Q_g: Gen103MVA.get_P_losses(P_g, Q_g, V_g).P_loss_tot

class TestModel1_derivatives(unittest.TestCase):  
    test_gen = Gen103MVA

    def test_general(self): 
        dP_loss_dQ = nd.Derivative(gen_loss_wrapper_Q_g, n=1)
        print(dP_loss_dQ())
    
if __name__ == "__main__": 
    unittest.main() 