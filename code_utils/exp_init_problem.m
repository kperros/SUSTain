function initP = exp_init_problem(X, params, init_method)
% initialization routines

if init_method==1 % random initialization of the "patient" factor and sampling through the input for the rest of the factors
    if (ismatrix(X))
        for i=1: length(params.modes)
            assert(isequal(params.modes{i}, 'discrete'));
            if i==2
                flag = false;
                while (~flag)
                    
                    pat_picked = randi(size(X, 1), params.R, 1);
                    initU{2} = full(X(pat_picked, :))';
                    
                    initU{2} = compute_init_sc_round(initU{2}, params.l_bounds(i), params.u_bounds(i));
                    flag = rank(initU{2}) == params.R;
                end
            else
                %initU{i} = randi([params.l_bounds(i) params.u_bounds(i)], size(X, i), params.R);
                initU{i} = zeros(size(X, i), params.R);
                for r = 1: params.R
                    p = randperm(size(X, i)); % similar to Chi and Kolda, 2011
                    nbig = round((1/params.R)*size(X, i));
                    initU{i}(p(1:nbig), r) = randi([params.l_bounds(i) params.u_bounds(i)], nbig, 1);
                end
                assert(rank(initU{i})==params.R); % assert that factor matrix is still full-rank
            end
        end
    else
        assert(ndims(X)==3);
        for i=1:3
            initU{i} = zeros(size(X, i), params.R);
            if (i==1)
                for r = 1: params.R
                    p = randperm(size(X, i)); % similar to Chi and Kolda, 2011
                    nbig = round((1/params.R)*size(X, i));
                    initU{i}(p(1:nbig), r) = randi([params.l_bounds(i) params.u_bounds(i)], nbig, 1);
                end
                assert(rank(initU{i})==params.R); % assert that factor matrix is still full-rank
            end
        end
        
        %initU{1} = randi([params.l_bounds(1) params.u_bounds(1)], size(X, 1), params.R);
        for j=1: params.R
            
            flag = true;
            while(flag)
                pat_picked = randi(size(X, 1)); % pick patients at random
                
                ind = X.subs(:, 1)==pat_picked;
                temp_mat = sparse(X.subs(ind, end - 1), X.subs(ind, end), X.vals(ind), size(X, 2), size(X, 3));
                
                % compute DX and other mode sums for this patient
                msum = sum(temp_mat, 2);
                [maxv, maxi] = max(msum);
                initU{3}(:, j) = full(temp_mat(maxi, :))';
                
                msum = sum(temp_mat, 1);
                [maxv, maxi] = max(msum);
                initU{2}(:, j) = full(temp_mat(:, maxi));
                
                if (rank(initU{2}(:, 1:j))< j || rank(initU{3}(:, 1:j))< j)
                    flag = true;
                else
                    flag = false;
                end
            end
        end
        
        for i=2:3
            initU{i} = compute_init_sc_round(initU{i}, params.l_bounds(i), params.u_bounds(i));
        end
        
    end
    
    lambda = ones(params.R, 1);
    initP = ktensor(lambda, initU);
    
else % pure random initialization
    assert(init_method==0);
    for i=1: length(params.modes)
        if (isequal(params.modes{i}, 'real'))
            assert(isinf(params.u_bounds(i)));
            initU{i} = rand(size(X, i), params.R);
        else
            assert(isequal(params.modes{i}, 'discrete'));
            assert(~isinf(params.u_bounds(i)));
            initU{i} = randi([params.l_bounds(i) params.u_bounds(i)], size(X, i), params.R);
        end
    end
    lambda = ones(params.R, 1);
    initP = ktensor(lambda, initU);
end