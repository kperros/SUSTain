function [h, mult]  = compute_init_sc_round(h, l_bound, up_bound)
% will compute scaled_rounding if exceeding upper bound
% used to initialize the factor matrices based on input values

assert(l_bound>=0);
for j=1: size(h, 2)
    maxval = max(h(:, j));
    if (maxval == 0)
        error('Caution, will divide by zero!');
    end
    if (maxval <= up_bound)
        continue;
    end
    
    mult(j) = up_bound/maxval;
    
    h(:, j) =  round(mult(j) .* h(:, j));
end

h(h < l_bound) = l_bound;
