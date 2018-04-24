clear all;
nn = 190;
load('dt_0.001000 damp_30.000000 N_190 k=80.000000 br=10.500000 cr=0.500000 var=1.mat');
max = 0;
for ii = 1:nn-1
    for jj = ii+1:nn
        delta = sqrt((x(ii,1)-x(jj,1))^2+(x(ii,2)-x(jj,2))^2+(x(ii,3)-x(jj,3))^2);
        if delta > max
            max = delta;
        end
    end
end
max