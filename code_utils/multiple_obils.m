function x = multiple_obils(B,y,l,u,Q,R,G)
% Adopted from MILES software to enable caching computations referring to
% the B matrix, which is common among all the systems solved.

% x = obils(B,y,l,u) produces the solution to the overdetermined box
% constrained integer least squares problem min_{l<=x<=u}||y-Bx||
%
% Inputs:
%    B - m-by-n real matrix with full column rank
%    y - m-dimensional real vector
%    l - n-dimensional integer vector, lower bound 
%    u - n-dimensional integer vector, upper bound, l < u    
%
% Output:
%    x - n-dimensional integer vector (in double precision) for
%        the optimal solution
%

% Subfunctions: obils_reduction, obils_search

% Main References: 
% [1] S. Breen and X.-W. Chang. Column Reordering for 
%     Box-Constrained Integer Least Squares Problems, 
%     Proceedings of IEEE GLOBECOM 2011, 6 pages.
% [2] X.-W. Chang and Q. Han. Solving Box-Constrained Integer Least 
%     Squares Problems, IEEE Transactions on Wireless Communications,  
%     7 (2008), pp. 277-287.

% Authors: Xiao-Wen Chang, www.cs.mcgill.ca/~chang 
%          Xiangyu Ren
% Copyright (c) 2015. Scientific Computing Lab, McGill University.
% April 2015. Last revision: December 2015

[m,n] = size(B);

if m ~= size(y,1) || size(y,2) ~= 1 || ...
      n ~= size(l,1) || size(l,2) ~= 1 || ...
      n ~= size(u,1) || size(u,2) ~= 1     % Input error
    error('Input arguments have a matrix dimension error!')
end

l = ceil(l); u = floor(u);  % To make it work with real bounds
for i = 1 : n
    if l(i) >= u(i)
        error('Invalid upper bound or lower bound');
    end
end

% Reduction
[R,y,l,u,p] = multiple_obils_reduction(B,y,l,u,Q,R,G);

% Search 
z = obils_search(R,y,l,u);

% Reorder z to obtain the optimal solution
x = zeros(n,1);
for i = 1 : n
    x(p(i)) = z(i); 
end
