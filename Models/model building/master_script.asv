close all;
clear all;

%%

%N_set = [ 50, 80, 100, 200 400];
N_set = [ 400, 450, 500 ];
dt_set = [ 0.001 ];
damp_set = [30];
k_set = [ 80 ];
%br_set = [5, 10, 13, 15, 17, 20 ];
br_set = [ 18 ];
core_set = [ 92 ];
globalvars.NN = 10;
globalvars.k_surf = 0.5;
%globalvars.damp = 2;
globalvars.rep_mult = 50;
globalvars.end_time = 500;
%globalvars.dt = 0.005;
%globalvars.coreR = 48;      % radius of PS sphere
%globalvars.eq_size = 13;    % Bead radius
var_set = [1 ];

%%
for j = 1:length(N_set)
    for i = 1:length(dt_set)
        for k = 1:length(damp_set)
            for l = 1:length(k_set)
                for m = 1:length(br_set )
                    for n = 1:length(core_set)
                        for p = 1:length(var_set)
                            globalvars.coreR = core_set(n);
                            globalvars.eq_size = br_set(m);
                            globalvars.dt = dt_set(i);
                            globalvars.k = k_set(l);
                            globalvars.damp = damp_set(k);
                            seed_centers_start = zeros( N_set(j), 3 );
                            ss = 0;
                            while ss < N_set(j)% randomly define initial bead position, < 300
                                trial =(rand( 1, 3 )-0.5)*2*300;
                                if euc_dist( trial, 0) < 300
                                    seed_centers_start( ss+1, : ) = trial;
                                    ss = ss+1;
                                end
                            end
                            tic   % count time
                            
                            
                            [beadList, stats] = BDAC_builder_driver_no_gfx( N_set(j), 'harm_trunc', seed_centers_start, globalvars );
                            %temp_struct = struct('pop_arr', pop, 'val', val_set(i), 'label', 'sigmoid' );
                            %struct_arr(i) = temp_struct;
                            toc
%                             
%                                         figure;
%                                         [xs,ys,zs] = sphere(20);
%                             
%                                         h = surf(48*xs,48*ys,48*zs);
%                                         set(h,'Facecolor',[1 0 0])
%                                         l = light;
%                                         lighting phong
%                                         hold on;
%                             
%                                         for  cc= 1:length(beadList)
%                             
%                                             h = surf(br_set(m)*(xs)+x(cc,1),br_set(m)*(ys)+x(cc,2),br_set(m)*(zs)+x(cc,3));
%                                             set(h,'Facecolor',[0 1 0])
%                                             plot3(x(cc,1), x(cc,2), x(cc,3), '.r' ); hold on;
%                                         end
%                             
                            h=figure;
                            plot( stats(:,1), stats(:,2), 'x' );
                            xlabel('Time');
                            ylabel('Average Velocity');
                            computed_title = sprintf('N=%d Dt=%f Damp=%f k=%f br=%f cr=%f var=%d', N_set(j), dt_set(i), damp_set(k), k_set(l), br_set(m), core_set(n), var_set(p) );
                            computed_png_file = sprintf('N=%d Dt=%f Damp=%f k=%f br=%f cr=%f var=%d.png', N_set(j), dt_set(i), damp_set(k), k_set(l), br_set(m), core_set(n), var_set(p) );
                            computed_mat_file = sprintf('N=%d Dt=%f Damp=%f k=%f br=%f cr=%f var=%d.fig', N_set(j), dt_set(i), damp_set(k), k_set(l), br_set(m), core_set(n), var_set(p) );
                            title(computed_title);
                            print( h, '-dpng', computed_png_file );
                            hgsave( computed_mat_file );
                            
                            for cct = 1:length(beadList)
                                x(cct, 1:3) = beadList(cct).x;
                                x_prime(cct, 1:3) = beadList(cct).x_prime;
                                v(cct, 1:3) = beadList(cct).v;
                                v_prime(cct, 1:3) = beadList(cct).v_prime;
                            end
                            
                            figure;
                            hist( euc_dist( x,0 ), 100);
                            savefile = sprintf('dt_%f damp_%f N_%d k=%f br=%f cr=%f var=%d.mat', globalvars.dt, globalvars.damp, N_set(j), k_set(l), br_set(m), core_set(n), var_set(p) );
                            save(savefile);
                            savefile = sprintf('dt_%f damp_%f N_%d k=%f br=%f cr=%f var=%d positions.mat', globalvars.dt, globalvars.damp, N_set(j), k_set(l), br_set(m), core_set(n), var_set(p) );
                            save( savefile, 'x', '-ascii', '-double' );
                            
                            clear dist_vec;
                            kk = 1;
                            for( ii = 1:length(x) )
                                for( jj = 1:length(x) )
                                    dist_vec(kk) = euc_dist( x(ii,:), x(jj,:));
                                    kk = kk + 1;
                                end
                            end
                            
                            figure;
                            hist( dist_vec, 1000 );
                        end
                    end
                end
                
            end
        end
    end
end
