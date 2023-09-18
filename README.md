# syngenlib

SynGenLib is a Python coded library for estimating the steady-state power losses in wound-rotor synchronous generators. The three models in this library differs in the following way: <br> 
Model 1 Uses [1] for its calculation model. Therefore, the Model 1 data class contains the nominal power losses for different generator components. <br> 
Model 2 uses resistance values in the stator and rotor at nominal operating point to estimate the rotor and stator losses. Otherwise, the core, stray, and other losses are calculated equivalently as Model 1. <br> 
Model 3 is equivalent to Model 2 with the main difference being that the resistances in stator and rotor is temperature dependend, and therefore has to be given as an input during the power loss calculations. 


### References 
[1] E. d. C. Bortoni, R. T. Siniscalchi, S. Vaschetto, M. A. Darmani, and A. Cavagnino, “Efficiency mapping and weighted average efficiency for large hydrogenerators,” IEEE Open J. Ind. Appl., vol. 2, pp. 11–20, 2021.
