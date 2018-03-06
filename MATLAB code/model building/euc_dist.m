function dist = euc_dist( x1, x2 )
    dist = sqrt( sum((x1 - x2).^2,2) ) + 1e-6;
end