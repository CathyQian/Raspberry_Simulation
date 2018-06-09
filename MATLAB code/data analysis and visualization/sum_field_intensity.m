
%this is to calculate the E4 of individual GNCs
%the only difference between individual GNC and dimer is that the dimer has
%"offset" parameter and therefore need to calculate dist1 and dist2
clear all
label1 = 'BDAC_N100_R48_r16_ps1_Air_785_3D_';
label2 = 'BDAC_N100_R48_r16_ps1_Air_860_3D_';

%dist_1 = 92e-9;
%dist_2 = 133.7e-9;
%dist_3 = 158.8e-9;
%dist_4 = 178.7223*10^-9;% for 184 nm core, 177 for 1 nm, 178 for 2 nm
% for 94 nm core, 110 for 1 nm, 111 for 2 nm
%dist_min = dist_3; %define E4 calculation region 
%dist_max = dist_4;


% for 94.5 nm core
dist_1 = 48e-9;
dist_2 = 90e-9;
dist_3 = 110e-9;% for 184 nm core, 177 for 1 nm, 178 for 2 nm
% for 94 nm core, 110 for 1 nm, 111 for 2 nm
dist_min = dist_2; %define E4 calculation region 
dist_max = dist_3;

load(sprintf('%sfieldx.mat', label1));
load(sprintf('%sfieldy.mat', label1));
load(sprintf('%sfieldz.mat', label1));
load(sprintf('%spos.mat', label1));
ex1=field_x;
ey1=field_y;
ez1=field_z;
e21=abs(ex1).^2+abs(ey1).^2+abs(ez1).^2;
sim_origin = [20.6017, 32.8531, 45.8893]*10^-9;% for 94 nm core
%sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;% for 184 nm core
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );
dist=sqrt(xmg.^2+ymg.^2+zmg.^2);

e21(dist>=dist_max)=0;
e21(dist<dist_min)=0;

load(sprintf('%sfieldx.mat', label2));
load(sprintf('%sfieldy.mat', label2));
load(sprintf('%sfieldz.mat', label2));

ex2=field_x;
ey2=field_y;
ez2=field_z;
e22=abs(ex2).^2+abs(ey2).^2+abs(ez2).^2;

e22(dist>=dist_max)=0;
e22(dist<dist_min)=0;

e2=e21.*e22;

sum(sum(sum(e2)))