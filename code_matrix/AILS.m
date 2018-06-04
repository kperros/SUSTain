function [U,V,lambda] = AILS(M, params)
U = params.P.U{1};
V = params.P.U{2}';
lambda = params.P.lambda;
if (params.lambda_on)
    assert(sum(lambda>=1)==length(lambda));
else
    assert(sum(lambda==1)==length(lambda));
end
R = params.R;

% Initialization
nM = norm(M,'fro')^2;
iter = 0;
e = [];
fit = [];

for i=1: length(params.modes)
    l_b(:, i) = params.l_bounds(i).*ones(R, 1);
    u_b(:, i) = params.u_bounds(i).*ones(R, 1);
end

while iter <= params.maxiter
    
    %% compute lambda
    if (params.lambda_on)
        new_lambda = compute_matrix_lambdas(M, U, V, R);
        lambda = new_lambda;
    end
    
    %% solve for U
    if (strcmp(params.modes{1}, 'discrete'))
        Unew = anls_obils(M', V'*diag(lambda), l_b(:, 1), u_b(:, 1))';
    else
        assert(sum(params.l_bounds==0)==length(params.l_bounds)); %make sure the lower bound for all modes is zero
        Unew = anls(M, V'*diag(lambda));
    end
    U = Unew;
    
    %% solve for V
    if (strcmp(params.modes{2}, 'discrete'))
        Vnew = anls_obils(M, U*diag(lambda), l_b(:, 2), u_b(:, 2));
    else
        assert(sum(params.l_bounds==0)==length(params.l_bounds)); %make sure the lower bound for all modes is zero
        Vnew = anls(M', U*diag(lambda))';
    end
    V = Vnew;
    
    %% convergence check and error monitoring
    e = [e compute_mf_normresidual(nM, M, U*diag(lambda), V)]; %e = [e sqrt( (nM-2*sum(sum(V.*A))+ sum(sum(B.*(V*V')))) )];
    fit = [fit (1 - (e(end) / sqrt(nM)))];
    if (iter>0)
        if (mod(iter, 1)==0)
            fprintf('AILS, iter: %d, current fit: %e\n', iter, fit(end));
        end
        if (abs(fit(end-1) - fit(end)) < params.stop_crit)
            disp('AILS fit for each iteration:');
            fit
            break;
        end
    end
    iter = iter + 1;
end

end

function new_lambda = compute_matrix_lambdas(M, U, V, R)
assert(ismatrix(M));
max_input_val = full(max(max(M)));
krp = khatrirao(V', U);
krp = rank_check_perturb(krp);
new_lambda = obils(krp, M(:), ones(R, 1), max_input_val.*ones(R, 1));
end

function U = anls(M, V) % non-negative alternating least squares
U = nnlsm_blockpivot(V, M')'; % nmf toolbox from Jingu Kim and Haesun Park
end

function [normresidual] = compute_mf_normresidual(normA_sq, A, W, H)
squared_residual = max(normA_sq - 2*trace(H*(A'*W)) + trace((W'*W)*(H*H')), 0);
normresidual = sqrt(squared_residual);
end
