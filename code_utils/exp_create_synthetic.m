function [X_true, X_true_ktensor, X] = exp_create_synthetic(dims, R_true, eta, u_bounds, l_bounds, lambda_up_bound)

for i=1: length(dims)
    assert(~isinf(u_bounds(i))); %\tau
    U{i} = zeros(dims(i), R_true);
    U{i} = randi([l_bounds(i) u_bounds(i)], dims(i), R_true);
end

lambda = randi([1 lambda_up_bound], R_true, 1);

X_true_ktensor = ktensor(lambda, U);
%X_true = full(X_true_ktensor); 
X_true = sptensor(full(X_true_ktensor));

if (ismatrix(X_true))
    X_true = double(X_true);
end

%X = X_true.*(sprand(dim1, dim2, sparsity)==0);%+ eta * (norm(X_true, 'fro')/norm(noise, 'fro')) * noise;
if (ismatrix(X_true))
    noise = rand(dims);
    X = X_true + (eta * (norm(X_true, 'fro')/norm(noise, 'fro')) * noise);
else
    assert(eta==0);
    X = X_true;
    %noise = tenrand(dims);
    %X = X_true + (eta * (norm(X_true)/norm(noise)) * noise);
end