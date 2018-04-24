function [x, x_prime, v, v_prime] = update_vect2(x, x_prime, v, v_prime, dt) % update our variables.

     xdif = euc_dist( x_prime, x );
     xdiftrig = xdif > 1;
     x( ~xdiftrig,: ) = x_prime( ~xdiftrig,: );
     x( xdiftrig,: ) = x(xdiftrig,:) + ( x_prime( xdiftrig,: ) - x(xdiftrig,:))./repmat(xdif(xdiftrig),[1,3]);
     

     vdif = euc_dist( v_prime, 0 );
     vdiftrig = vdif > 1;     
     v( ~vdiftrig,:) = v_prime( ~vdiftrig,: );     
     v( vdiftrig, : ) = v_prime( vdiftrig,: )./repmat(vdif(vdiftrig),1,3);

    
end
