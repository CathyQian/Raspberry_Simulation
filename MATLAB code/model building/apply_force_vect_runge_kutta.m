function[x, x_prime, v, v_prime] = apply_force_vect_runge_kutta( x, x_prime, v, v_prime, force, dt, globalvars )

stop
I'm broken

x_prime = x + v*dt + 0.5*force( x,  v,  globalvars)*dt*dt;

v_prime = v + dt*force( x, v,  globalvars);

    
end      