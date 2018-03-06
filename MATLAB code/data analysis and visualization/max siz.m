clear all;
load('dt_0.001000 damp_30.000000 N_800 k=80.000000 br=14.000000 cr=92.000000 var=1 positions.mat');
x = [x(:,1), x(:,2), x(:,3)];
max_f=100;
for i = 1:800
    dist=sqrt(x(i,1).^2+x(i,2).^2+x(i,3).^2);
        if  dist>max_f
            max_f=dist;
        end
end
