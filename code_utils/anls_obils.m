function V = anls_obils(M, U, l_b, u_b)
%% caching computations common for all systems solved
U = rank_check_perturb(U);
% Transform B and y by using the QR factorization
[Q,R] = qr(U,0);
% Inverse transpose of R
G = inv(R)';

%% solve with box-constrained integer least squares
for col=1: size(M, 2)
    V(:, col) = multiple_obils(U, M(:, col), l_b, u_b, Q, R, G);
    %assert(isequal(V(:, col), obils(U, M(:, col), l_b, u_b)));
end
% for col=1: size(M, 2)
%     assert(isequal(V(:, col), obils(U, M(:, col), l_b, u_b)));
% end
end