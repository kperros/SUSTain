function [V, lambda] = SUSTain_Update_Factor(V,coeff,mkrp,mode,lambda,params)
assert(strcmp(params.modes{mode}, 'discrete')); % assert that the mode has discrete constraint
assert(params.lambda_on == 1);

if (params.max_inner_iter>1)
    V_scaled_prev = V*diag(lambda); % accelerated version (the #inner iterations is decided based on the difference of the learnt matrix)
    init_diff_V = [];
end

r = size(V, 2);
for myiter=1: params.max_inner_iter
    for k = 1 : r        
        core = V*(lambda.*coeff(:, k)); % most expensive part of updating k-th column      
        core_k = V(:, k) .* (lambda(k)*coeff(k,k)); % will be subtracted to update based on new lambda(k)
        
        %% update scalar weight of k-th component
        delta_lambda_k = ( V(:, k)' * ( mkrp(:, k) - core ) )/( coeff(k,k) * norm(V(:, k))^2 );
        lambda(k) = max(1, round(lambda(k) + delta_lambda_k));
        
        core = core - core_k + ( V(:, k) .* (lambda(k)*coeff(k,k)) ); % adjust core with new lambda(k)
        
        %% update k-th column of factor matrix
        Vk = V(:,k) + ( (mkrp(:, k) - core)/(coeff(k,k)*lambda(k)) );
        
        %% project to the integer constrained set
        Vk = round(Vk);
        if (~isinf(params.l_bounds(mode)))
            Vk(Vk < params.l_bounds(mode)) = params.l_bounds(mode);
        end
        if (~isinf(params.u_bounds(mode)))
            Vk(Vk > params.u_bounds(mode)) = params.u_bounds(mode);
        end
        V(:,k) = Vk;
        
        %% avoid zero-lock
        if V(:,k) == 0
            assert(params.u_bounds(mode) >= 1 && params.l_bounds(mode)>=0);            
            V(randi([1 size(V, 1)], 1, 1), k) = 1; % smallest perturbation possible
        end
    end
    
    if (params.max_inner_iter>1)
        if (myiter==1)
            init_diff_V = norm(V*diag(lambda) - V_scaled_prev, 'fro');
        else
            if (norm(V*diag(lambda) - V_scaled_prev, 'fro') <= params.alpha * init_diff_V )
                break;
            end
        end
        V_scaled_prev = V*diag(lambda);
    end
end
