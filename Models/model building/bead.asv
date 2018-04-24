classdef bead < handle 
   properties  % Dynamic Properties
      x = [0,0,0]; % store current step info
      v = [0,0,0]; 

      id = -1;
      
      x_eq = 1; % Bead radius
      
      x_prime = [0,0,0]; % store 1 step ahead info
      v_prime = [0,0,0];
      
      time_since_birth = 0;
      time_last_divide = 0;
      
      gen = 0;
   end
   
   methods
      function cc = bead(x_init, v_init, id_init, beadR) %construcor
         cc.x = x_init;
         cc.v = v_init;
         cc.id = id_init;
         cc.x_eq = beadR;
      end
      
      function update(cc, dt) % update our variables.
          
          % Max movement rule
          if( euc_dist( cc.x_prime, cc.x ) > 1 )
              cc.x = cc.x + (cc.x_prime - cc.x )/euc_dist( cc.x_prime, cc.x);
          else
              cc.x = cc.x_prime;
          end
          
          % Max velocity rule
          if( euc_dist(cc.v_prime,0) > 1 )
              cc.v = cc.v_prime*1/euc_dist(cc.v_prime,0);
          else
              cc.v = cc.v_prime;
          end
      end
      

      
      % Force is a function handle to something that takes an object, and
      % an object list. More or less we can tune force to do what we want
      % instead of flubbing around with the crypt cell function.
      % Note that this is roughly the eulor method. It most certainly isn't
      % stable. But if we keep the force highly damped, I don't think we
      % really care. We don't need the exact phase space progression, just
      % a plausable one, which this will find.    
      function apply_force( cc, cc_ind, force, global_cells, dt, globalvars )
                
          
            cc.v_prime = cc.v + dt*force(cc, cc.x, cc_ind, global_cells, globalvars );
            
            if euc_dist( cc.v_prime, 0 ) > 1
                cc.v_prime = (cc.v_prime )/euc_dist( cc.v_prime, 0 );
            end                
            
            cc.x_prime = cc.x + (cc.v_prime + cc.v)/2*dt;                      
        
            if cc.x_prime(1) > globalvars.bounds(1,2) || cc.x_prime(1) < globalvars.bounds(1,1)
                cc.x_prime(1) = cc.x(1);
            end
            if cc.x_prime(2) > globalvars.bounds(2,2) || cc.x_prime(2) < globalvars.bounds(2,1)
                cc.x_prime(2) = cc.x(2);
            end
            if cc.x_prime(3) > globalvars.bounds(3,2) || cc.x_prime(3) < globalvars.bounds(3,1)
                cc.x_prime(3) = cc.x(3);
            end     
     

      end      
   end      
end