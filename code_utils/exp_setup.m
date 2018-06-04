function params = exp_setup(u_bounds, l_bounds, maxiter, stop_crit)

if (length(u_bounds)==3)
    params.modes = {'discrete', 'discrete', 'discrete'};
else
    assert(length(u_bounds)==2);
    params.modes = {'discrete', 'discrete'};
end

params.lambda_on = 1;
params.u_bounds = u_bounds;
params.l_bounds = l_bounds;
params.maxiter = maxiter;
params.stop_crit = stop_crit;
params.max_inner_iter = 1;
params.alpha = 0.01; % threshold for breaking inner iterations if params.max_inner_iter > 1