n = 100;% number of beads
load('dt_0.001000 damp_30.000000 N_800 k=80.000000 br=11.000000 cr=92.000000 var=1.mat');
i = 1;
x1 = [x(:,1), x(:,2), x(:,3)];
xt = zeros(400,3);
%xt = x1;
%xt(i) = [x(j,1), x(j,2), x(j,3)]; and start with i=1
for j = 1:800
    if x1(j,1)>0
            % xt[x(i,1), x(i,2), x(i,3)] = [x(j,1), x(j,2), x(j,3)];
            xt(i,1) = x(j,1); 
            xt(i,2) = x(j,2);
            xt(i,3) = x(j,3);
            i = i+1;
    end
end

savefile = sprintf('half_medium_GNC.mat');
save( savefile, 'xt', '-ascii', '-double' ); % for Lumerical FDTD loading bead position
savefile = sprintf('half_medium_GNC_copy.mat');
save( savefile, 'xt', '-mat', '-double' ); % for Matlab data processing