clear all; close all;
addpath('../code_utils/');

%% params
R = 10;
u_bounds = [5 5 5];
l_bounds = [0 0 0];
assert(sum(u_bounds>0)==length(u_bounds));

maxiter = 1e6;
stop_crit = 1e-4; % successive difference of the objective so that we break iterations

eta_noise = 0;
lambda_up_bound = 50; % upper bound for lambda generation
true_rank = 50;
synth_sz = 50;
[X_true, X_true_ktensor, X] = exp_create_synthetic([synth_sz, synth_sz, synth_sz], true_rank, eta_noise, u_bounds,l_bounds, lambda_up_bound);

params = exp_setup(u_bounds, l_bounds, maxiter, stop_crit);
params.R = 10;

%% methods
% initial point
params.P = exp_init_problem(X, params, 0);

% SUSTain_T
P = SUSTain_T(X, params);