The current repository provides the code accompanying our KDD 2018 paper "SUSTain: Scalable Unsupervised Scoring for Tensors and its Application to Phenotyping", by Ioakeim Perros, Evangelos E. Papalexakis, Haesun Park, Richard Vuduc, Xiaowei Yan, Christopher deFilippi, Walter F. Stewart, Jimeng Sun.

For each of the matrix and tensor versions of SUSTain, there exist two separate folders named "code_matrix" and "code_tensor" respectively. Also, there is a folder named "code_utils" containing utility functions useful for both matrix and general tensor versions.
Within the folder "code_matrix", the starting point is the file "exp_matrix.m".
Within the folder "code_tensor", the starting point is the file "exp_tensor.m".

We have tested the code on MatlabR2017b. The prerequisite packages to run it are:
- Tensor Toolbox Version 2.6 (can be downloaded from: http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html)
- MILES library for solving box-constrained Integer Least Squares problems (can be downloaded from: http://cs.mcgill.ca/~chang/MILES_routine3.php)

List of files included: 
code_matrix: exp_matrix.m, SUSTain_M.m, AILS.m

code_tensor: exp_tensor.m, SUSTain_T.m

code_utils: SUSTain_Update_Factor.m, exp_create_synthetic.m, exp_init_problem.m, exp_setup.m, anls_obils.m, rank_check_perturb.m, multiple_obils.m, multiple_obils_reduction.m, compute_init_sc_round.m