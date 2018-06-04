function P = SUSTain_T(X, params)
% SUSTain_T procedure

printitn = 1; % after how many iters we will print fit information
N = ndims(X);
normX = norm(X);
U = params.P.U;
lambda = params.P.lambda;
R = params.R;
dimorder = 1:N;
fit = 0;

%% init matrix of coefficients
UtU = zeros(R,R,N);
for n = 1:N
    if ~isempty(U{n})
        UtU(:,:,n) = U{n}'*U{n};
    end
end

%% Main loop
for iter=1: params.maxiter
    fitold = fit;
    %% Iterate over all N modes of the tensor
    for n = dimorder(1:end)
        
        %% Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
        Unew = mttkrp(X,U,n); % tensor toolbox (Kolda et al.)
        
        %% Compute the matrix of coefficients for linear system
        Y = prod(UtU(:,:,[1:n-1 n+1:N]),3);
        % Unew = Unew / Y;
        
        %% SUSTain_Update_Factor updating n-th factor matrix and lambdas
        [Unew, lambda] = SUSTain_Update_Factor(U{n},Y,Unew,n,lambda,params);
        
        U{n} = Unew;
        UtU(:,:,n) = U{n}'*U{n};
    end
    
    P = ktensor(lambda,U);
    if normX == 0
        fit = norm(P)^2 - 2 * innerprod(X,P);
    else
        normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
        fit = 1 - (normresidual / normX); %fraction explained by model
    end
    fitchange = abs(fitold - fit);
    
    % Check for convergence
    if (iter > 1) && (fitchange < params.stop_crit)
        flag = 0;
    else
        flag = 1;
    end
    
    if (mod(iter,printitn)==0) || ((printitn>0) && (flag==0))
        fprintf('SUSTain_T: Iter %2d: f = %e f-delta = %7.1e\n', iter, fit, fitchange);
    end
    
    % Check for convergence
    if (flag == 0)
        fprintf('SUSTain_T done.\n\n');
        break;
    end
end
