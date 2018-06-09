%% integral all area of smooth nanoshell dimers
clear all
label1 = 'Smooth NS_R92_R169_Air_3D_785_';
label2 = 'Smooth NS_R92_R169_Air_3D_860_';

load(sprintf('%sfieldx.mat', label1));
load(sprintf('%sfieldy.mat', label1));
load(sprintf('%sfieldz.mat', label1));
load(sprintf('%spos.mat', label1));
ex1=field_x;
ey1=field_y;
ez1=field_z;
e21=abs(ex1).^2+abs(ey1).^2+abs(ez1).^2;

load(sprintf('%sfieldx.mat', label2));
load(sprintf('%sfieldy.mat', label2));
load(sprintf('%sfieldz.mat', label2));

ex2=field_x;
ey2=field_y;
ez2=field_z;
e22=abs(ex2).^2+abs(ey2).^2+abs(ez2).^2;

e2 = e21.*e22;

sum(sum(sum(e2)))
%% integral over 1 nm or 2 nm surface of nanoshell dimers
clear all
label1 = 'Smooth NS_R92_R169_dimer_sep_1nm_Air_3D_785_';
label2 = 'Smooth NS_R92_R169_dimer_sep_1nm_Air_3D_860_';

'Dont forget to change offset1!!!'

offset = 169.5*10^-9;% 169.5 for 1 nm separation, 171.5 nm for 5 nm separation, 419 for 500 nm separation
inter_scale = 92*10^-9;
max_scale = 171*10^-9;%170 for 1 nm, 171 for 2 nm

load(sprintf('%sfieldx.mat', label1));
load(sprintf('%sfieldy.mat', label1));
load(sprintf('%sfieldz.mat', label1));
load(sprintf('%spos.mat', label1));
ex1=field_x;
ey1=field_y;
ez1=field_z;
e21=abs(ex1).^2+abs(ey1).^2+abs(ez1).^2;

sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;% smooth nanoshells, same as 184 nm core result
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );
dist1=sqrt((xmg+offset).^2+ymg.^2+zmg.^2);
dist2=sqrt((xmg-offset).^2+ymg.^2+zmg.^2);
dist0=sqrt(xmg.^2+ymg.^2+zmg.^2);

e2(~(((dist1<=max_scale)&(dist1>inter_scale))|((dist2<=max_scale)&(dist2>inter_scale))))=0;

e21(~(((dist1<=max_scale)&(dist1>inter_scale))|((dist2<=max_scale)&(dist2>inter_scale))))=0;


load(sprintf('%sfieldx.mat', label2));
load(sprintf('%sfieldy.mat', label2));
load(sprintf('%sfieldz.mat', label2));

ex2=field_x;
ey2=field_y;
ez2=field_z;
e22=abs(ex2).^2+abs(ey2).^2+abs(ez2).^2;

e22(~(((dist1<=max_scale)&(dist1>inter_scale))|((dist2<=max_scale)&(dist2>inter_scale))))=0;

e2 = e21.*e22;

sum(sum(sum(e2)))
%% integral over 1 nm or 2 nm of smooth nanoshell surface
clear all
label1 = 'Smooth NS_R92_R169_Air_3D_785_';
label2 = 'Smooth NS_R92_R169_Air_3D_860_';

'Dont forget to change offset2!!!'

dist_min = 92*10^-9;
dist_max = 171*10^-9;%170 for 1 nm, 171 for 2 nm

load(sprintf('%sfieldx.mat', label1));
load(sprintf('%sfieldy.mat', label1));
load(sprintf('%sfieldz.mat', label1));
load(sprintf('%spos.mat', label1));
ex1=field_x;
ey1=field_y;
ez1=field_z;
e21=abs(ex1).^2+abs(ey1).^2+abs(ez1).^2;

%sim_origin = [20.6017, 32.8531, 45.8893]*10^-9;% for 94 nm core
sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;% for 184 nm core and smooth nanoshell
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );
dist0=sqrt(xmg.^2+ymg.^2+zmg.^2);

e2(~(((dist0<=dist_max)&(dist0>dist_min))))=0;
e21(~(((dist0<=dist_max)&(dist0>dist_min))))=0;

load(sprintf('%sfieldx.mat', label2));
load(sprintf('%sfieldy.mat', label2));
load(sprintf('%sfieldz.mat', label2));
ex2=field_x;
ey2=field_y;
ez2=field_z;
e22=abs(ex2).^2+abs(ey2).^2+abs(ez2).^2;

e21(~(((dist0<=dist_max)&(dist0>dist_min))))=0;

e2 = e21.*e22;

sum(sum(sum(e2)))
%%

