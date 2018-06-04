function [U,V,lambda] = SUSTain_M(X, params)
% SUSTain_M procedure

U = params.P.U{1};
V = params.P.U{2};
lambda = params.P.lambda;
if (params.lambda_on)
    assert(sum(lambda>=1)==length(lambda));
else
    assert(sum(lambda==1)==length(lambda));
end
nX = norm(X,'fro')^2;
iter = 0;
e = [];
fit = [];

%% Main loop
while iter <= params.maxiter
    %% Update of U and lambdas
    A = X*V; 
    B = V'*V;
    [U, lambda] = SUSTain_Update_Factor(U,B,A,1,lambda,params);
    
    %% Update of V and lambdas
    A = X'*U; 
    B = U'*U;
    [V, lambda] = SUSTain_Update_Factor(V,B,A,2,lambda,params);
    
    %% convergence check and error monitoring
    %fprintf('Iter %d: after updating V, residual norm: %e \n',iter,compute_mf_normresidual(nM, X, U, V));
    % e = [e compute_mf_normresidual(nX, X, U*diag(lambda), V')];
    V_ = V*diag(lambda);
    e = [e sqrt( max(nX- 2*sum(sum(V_.*A)) + sum(sum(B.*(V_'*V_))), 0) ) ];    
    fit = [fit (1 - (e(end) / sqrt(nX)))];
    if (iter>0)
        if (abs(fit(end-1) - fit(end)) < params.stop_crit)
            disp('SUSTain_M fit for each iteration:');
            fit
            break;
        end
    end
    iter = iter + 1;
end
