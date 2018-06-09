function of = nn_force_vect2( x,  v,  globalvars)
    k = globalvars.k;
    k_surf = globalvars.k_surf;        

    damp = globalvars.damp;
    rep_mult = globalvars.rep_mult;
    
    of = zeros( size(x));              
    
    
    NN = globalvars.NN;
    
    vectDif = zeros( length(x), 3, NN );
    nDist3d= zeros( length(x), 3, NN );
    nDistDifs3d = zeros( length(x), 3, NN );
    
    
    for i = 1:NN
        vectDif(:,:,i) = reshape( x - x(globalvars.neighborIds(:,i),:), length(x), 3, 1);
        nDist3d(:,:,i) = repmat( euc_dist(x, x(globalvars.neighborIds(:,i),:)), [1, 3, 1]);
        nDistDifs3d(:,:,i) = repmat( euc_dist(x, x(globalvars.neighborIds(:,i),:))-2*globalvars.eq_size, [1, 3, 1]);
    end
    
    nDistWeights = nDist3d < 2*globalvars.eq_size;
    
    %size( vectDif)
    %size( nDist3d )
    %size( nDistDifs3d )
    
    force_expanded = k*vectDif./nDist3d.*nDistDifs3d.*nDistWeights;
    force_expanded(isnan(force_expanded)) = 0;
    force_repel = sum( force_expanded, 3 );
    %size(force_repel)
    
    %%    % Harmonic attraction to surface
    surf_pos = x./repmat(euc_dist( x, 0 ),1,3)*(globalvars.coreR + globalvars.eq_size);
    %surf_dist = euc_dist( x, surf_pos );
    surf_mult_bin = (euc_dist( x, 0 ) >  euc_dist( surf_pos, 0 )) + (euc_dist( x, 0 ) <=  euc_dist( surf_pos, 0 ))*rep_mult;
    surf_mult = repmat( surf_mult_bin, 1, 3 );
    surf_force = k_surf*(x - surf_pos).*surf_mult;
    
    of = of - surf_force;
    of = of - force_repel;
    of = of - damp*v;
    
end