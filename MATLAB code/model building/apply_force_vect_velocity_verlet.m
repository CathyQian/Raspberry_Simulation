function[x, x_prime, v, v_prime, at_prime] = apply_force_vect_velocity_verlet( x,  v,  force, dt, globalvars, at, carry )

% VV doesn't actually handle damping properly, but we're going to ignore
% that. 
if( carry < 1 )
    at = force( x,  v,  globalvars);
end

x_prime = x + v*dt + 0.5*at*dt*dt;

at_prime = force( x_prime, v, globalvars );

v_prime = v + dt*0.5*(at+at_prime);
    
end      