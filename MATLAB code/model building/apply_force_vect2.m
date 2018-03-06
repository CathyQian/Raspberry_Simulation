function[x, x_prime, v, v_prime] = apply_force_vect2( x, x_prime, v, v_prime, force, dt, globalvars )

v_prime = v + dt*force( x,  v,  globalvars);

vdif = euc_dist( v_prime, 0 );
vdiftrig = vdif > 1;
v_prime( vdiftrig, : ) = v_prime( vdiftrig,: )./repmat(vdif(vdiftrig),1,3);

x_prime = x + (v_prime + v)/2*dt;

    
end      