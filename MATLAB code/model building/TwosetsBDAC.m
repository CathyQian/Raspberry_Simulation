d = 1;
n = 100;
load('');
x1 = [x(:,1) - d/2, x(:,2), x(:,3)];
x2 = [x(:,1) + d/2, x(:,2), x(:,3)];
xt = [x1;x2];
min = 100000;
for i = 1:n
    for j = 1:n
        d = sqrt((x1(i,1)-x2(j,1))^2+(x1(i,2)-x2(j,2))^2+(x1(i,3)-x2(j,3))^2);
        if d < min
            min = d;
        end
    end
end
min
%savefile = 'Two_BDAC_offset.mat';
%save( savefile, 'xt', '-ascii', '-double' );