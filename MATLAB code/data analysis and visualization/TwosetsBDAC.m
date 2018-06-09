dist = 338;% left or right shift
n = 800;% number of beads
insertname=int2str(dist);
load('dt_0.001000 damp_30.000000 N_800 k=80.000000 br=14.000000 cr=92.000000 var=1.mat');
x1 = [x(:,1) - dist/2, x(:,2), x(:,3)];
x2 = [x(:,1) + dist/2, x(:,2), x(:,3)];
min_1 = [0,0,0];
min_2 = [0,0,0];
xt = [x1;x2];
min = dist;
label = 0;
for i = 1:800 % ignore i=0 which is PS beads
    for j = 1:800
        if x1(i,1)>x2(j,1)
            label = 1;% check whether the two clusters overlap
        end
        d = sqrt((x1(i,1)-x2(j,1))^2+(x1(i,2)-x2(j,2))^2+(x1(i,3)-x2(j,3))^2);
        if d < min
            min = d;
            min_1 = [x1(i,1), x1(i,2), x1(i,3)];
            min_2 = [x2(j,1), x2(j,2), x2(j,3)];
        end
    end
end
min = min - 26;

savefile = sprintf('Two_BDAC_offset_dist_%s_184core_GNC_sep_9nm.mat',insertname);
save( savefile, 'xt', '-ascii', '-double' ); % for Lumerical FDTD loading bead position
savefile = sprintf('Two_BDAC_offset_dist_%s_184core_GNC_sep_9nm_copy.mat',insertname);
save( savefile, 'xt', '-mat', '-double' ); % for Matlab data processing
