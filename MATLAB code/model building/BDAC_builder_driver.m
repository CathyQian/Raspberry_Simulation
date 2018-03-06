function [pop, stat] = BDAC_builder_driver( val, label, seed_centers )

% setup simulation parameters
end_time = 70;    % end time
dt = 0.001;                % time step
steps = end_time / dt;  % we'll use linear steps
NN = 7;                 % Number of neighbors to consider
globalvars.NN=NN;

globalvars.centers = seed_centers;

globalvars.val = val;

%units of distance is 1 nm
globalvars.coreR = 48;      % radius of PS sphere
globalvars.eq_size = 13;    % Bead radius
[beads,trash] = size(globalvars.centers);
globalvars.bounds = [-10000,10000;-10000,10000; -10000,10000 ]; % 20 micron size
   
for i = 1:beads
    beadList(i) = bead(globalvars.centers(i,:), [0,0,0], i, globalvars.eq_size); %#ok<AGROW>
end

id = i+1;



%record = zeros( steps, 3*beads);

file3d = sprintf('beads3d_%s_%d.avi', label, val );
file3d %#ok<NOPRT>
aviobj3d=avifile(file3d); %creates AVI file

steps_since_recalc = 1e6;
k = 1;
    
for cct = 1:beads
    x(cct, 1:3) = beadList(cct).x;
    x_prime(cct, 1:3) = beadList(cct).x_prime;
    v(cct, 1:3) = beadList(cct).v;
    v_prime(cct, 1:3) = beadList(cct).v_prime;
end

for i = 1:steps
        
    
    if (beads > 1 && steps_since_recalc > 20 )
        % Work out the NN. The KNN code used is borrowed, see file for credits.
        if beads < NN
            use_nn = length(beadList);
        else
            use_nn = NN;
        end
        
        dataMatrix = zeros( beads, 3 );
        

        dataMatrix = x;

        
        queryMatrix = dataMatrix;        
        
        [globalvars.neighborIds globalvars.neighborDistances] = kNearestNeighbors(dataMatrix, queryMatrix, use_nn);
        %globalvars.neighborIds
        %stop
        steps_since_recalc = 0;                
    else
        steps_since_recalc = steps_since_recalc+1;
    end

                     
        
    [x, x_prime, v, v_prime] = apply_force_vect2( x, x_prime, v, v_prime, @nn_force_vect2, dt, globalvars );        
    [x, x_prime, v, v_prime] = update_vect2(x, x_prime, v, v_prime, dt);
    
    
    if mod(i,200)==0               
        fprintf('%d of %d\n', i, steps );

                
        hf= figure('visible','off');
        hax=axes;     
                
        vel = sum(euc_dist(v,0));
                    
        plot3(x(:,1), x(:,2), x(:,3), '.k', 'parent', hax ); hold on;        

        set(gca,'Xlim',[-500 500]); 
        set(gca,'Ylim',[-500 500]);    
        set(gca, 'Zlim', [ -500 500]);
        aviobj3d = addframe(aviobj3d, hf);
        close(hf); %closes the handle to invisible figure close all;   

        stat_gen(k,1:2) = [ i*dt, vel ];
        k = k + 1;                
    end
    
    
    
end

aviobj3d=close(aviobj3d); %closes the AVI file

%tic
%M = make_movie( 'test.avi', record );
%toc

% Not that this makes a lot of sense at the moment, but here is the
% feedback.

for cct = 1:beads
    beadList(cct).x = x(cct, 1:3);
    beadList(cct).x_prime = x_prime(cct, 1:3);
    beadList(cct).v = v(cct, 1:3);
    beadList(cct).v_prime = v_prime(cct, 1:3);
end

pop = beadList;
stat = stat_gen/beads;

% figure;
% plot( stat_gen(:,1), stat_gen(:,2), 'xr' );
% figure;
% plot( stat_total_cells(:,1), stat_total_cells(:,2), 'ob' );

clear beadList
