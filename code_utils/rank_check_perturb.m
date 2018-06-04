function Unew =  rank_check_perturb(U)
if rank(U) < min(size(U))
    
    Unew = U;
    pwr = -12;
    iter = 1;
    while (rank(Unew) ~= min(size(Unew)) && iter<10)
        Unew = U + power(10, pwr).*eye(size(U));
        fprintf('Perturbed the rank of the integer LS matrix. Just used: %e\n',power(10, pwr));
        
        pwr = pwr + 2;
        iter = iter + 1;
    end
    assert(rank(Unew) == min(size(Unew)));
    disp('End of perturbing.');
    %error('Matrix does not have full column rank, use ubils')
else
    Unew = U;
end