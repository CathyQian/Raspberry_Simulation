function[x, x_prime, v, v_prime] = apply_force_vect_2nd_order_euler( x,  v,  appliedForce, dt, globalvars )

f = appliedForce( x,  v,  globalvars);

v_prime = v + dt*f;
x_prime = x + dt*v_prime + 0.5*f*dt*dt;
    
end      