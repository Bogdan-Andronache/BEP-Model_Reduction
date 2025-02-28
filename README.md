# PoissonModelReduction
"# BEP-Model_Reduction" 
"# BEP-Model_Reduction" 
The 3 folders called NonLinearPoisson followed by the parameter space contain the newest code. Each of them contains scripts which I used to generate parameters and snapshots. There is also a folder with results. There are 3 main model scripts which make use of the EZyRB library: model_rank.py (uses RBF to plot error against the rank used for SVD), model_sampling.py (plots the maximum error found in the mesh for each testing sanpshot in the right grid position) and model_singular_values_modes.py (plots the first 30 singular_values and their respective modes). There is also a script called model_ANN.py which does the same as model_rank but using ANN.
