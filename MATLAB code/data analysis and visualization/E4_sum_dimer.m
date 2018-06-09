%% BDAC nanoshell dimer (separate E at 785 nm and 860 nm)
clear all
label1 = 'BDAC_N800_R92_r14_ps1_dimer_sep_5nm_Air_3D_785_p_';
label2 = 'BDAC_N800_R92_r14_ps1_dimer_sep_5nm_Air_3D_860_p_';

'Dont forget to change offset!!!'

offset = ((329+5)*10^-9)/2;

% for 184 nm core
dist_1 = 92e-9;
dist_4 = 178.7223*10^-9;% for 184 nm core, 177 for 1 nm, 178 for 2 nm
% for 94 nm core, 110 for 1 nm, 111 for 2 nm
dist_min = dist_1; %define E4 calculation region 
dist_max = dist_4;

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

sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;% for 184 nm core
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );
dist1=sqrt((xmg+offset).^2+ymg.^2+zmg.^2);
dist2=sqrt((xmg-offset).^2+ymg.^2+zmg.^2);

e2(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;
e21(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;
e22(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;

e2 = e21.*e22;

sum(sum(sum(e2)))

%% BDAC nanoshell dimer (E at 785 nm and 860 nm in the same file)
clear all
label1 = 'BDAC_N800_R92_r14_ps1_dimer_sep_2nm_Air_3D_785_860_';

'Dont forget to change offset!!!'

offset = ((329+2)*10^-9)/2;
dist_min = 92*10^-9; %define E2 calculation region 
dist_max = 178.7223*10^-9;% for 184 nm core, 177 for 1 nm, 178 for 2 nm

load(sprintf('%sfieldx.mat', label1));
load(sprintf('%sfieldy.mat', label1));
load(sprintf('%sfieldz.mat', label1));
load(sprintf('%spos.mat', label1));

ex1=field_x(1:592,1:301,1:301,2);% 785 nm data
ey1=field_y(1:592,1:301,1:301,2);% 785 nm data
ez1=field_z(1:592,1:301,1:301,2);% 785 nm data
%wavelength_1 = 3e8/field_f(2)
e21=abs(ex1).^2+abs(ey1).^2+abs(ez1).^2;

ex2=field_x(1:592,1:301,1:301,1);% 860 nm data
ey2=field_y(1:592,1:301,1:301,1);% 860 nm data
ez2=field_z(1:592,1:301,1:301,1);% 860 nm data
%wavelength_2 = 3e8/field_f(1)
e22=abs(ex2).^2+abs(ey2).^2+abs(ez2).^2;

sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;% for 184 nm core
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );
dist1=sqrt((xmg+offset).^2+ymg.^2+zmg.^2);
dist2=sqrt((xmg-offset).^2+ymg.^2+zmg.^2);

e2(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;
e21(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;
e22(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;

e2 = e21.*e22;

sum(sum(sum(e2)))
%% Spherical Au NP dimer
clear all
label1 = 'AuSphere_R169_dimer_sep50nm_Air_3D_785_860';

'Dont forget to change offset!!!'

offset = ((338+50)*10^-9)/2;
dist_min = 92e-9; %define E2 calculation region 
dist_max = 171e-9; % 170 for 1 nm over surface, 171 for 2 nm over surface 

load(sprintf('%sfieldx.mat', label1));
load(sprintf('%sfieldy.mat', label1));
load(sprintf('%sfieldz.mat', label1));
load(sprintf('%spos.mat', label1));

ex1=field_x(1:396,1:181,1:181,2);% 785 nm data
ey1=field_y(1:396,1:181,1:181,2);% 785 nm data
ez1=field_z(1:396,1:181,1:181,2);% 785 nm data
%wavelength_1 = 3e8/field_f(2)
e21=abs(ex1).^2+abs(ey1).^2+abs(ez1).^2;

ex2=field_x(1:396,1:181,1:181,1);% 860 nm data
ey2=field_x(1:396,1:181,1:181,1);% 860 nm data
ez2=field_x(1:396,1:181,1:181,1);% 860 nm data
%wavelength_2 = 3e8/field_f(1)
e22=abs(ex2).^2+abs(ey2).^2+abs(ez2).^2;

sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;% for nanosphere dimers
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );
dist1=sqrt((xmg+offset).^2+ymg.^2+zmg.^2);
dist2=sqrt((xmg-offset).^2+ymg.^2+zmg.^2);

e2(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;
e21(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;
e22(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;

e2 = e21.*e22;

sum(sum(sum(e2)))

%% Smooth Au nanoshell dimer
clear all
label1 = 'Smooth NS_R92_R169_dimer_sep50nm_Air_3D_785_860';

'Dont forget to change offset!!!'

offset = ((338+50)*10^-9)/2;
dist_min = 92e-9; %define E2 calculation region 
dist_max = 171e-9; % 170 for 1 nm over surface, 171 for 2 nm over surface 

load(sprintf('%sfieldx.mat', label1));
load(sprintf('%sfieldy.mat', label1));
load(sprintf('%sfieldz.mat', label1));
load(sprintf('%spos.mat', label1));

ex1=field_x(1:401,1:181,1:181,2);% 785 nm data
ey1=field_y(1:401,1:181,1:181,2);% 785 nm data
ez1=field_z(1:401,1:181,1:181,2);% 785 nm data
%wavelength_1 = 3e8/field_f(2)
e21=abs(ex1).^2+abs(ey1).^2+abs(ez1).^2;

ex2=field_x(1:401,1:181,1:181,1);% 860 nm data
ey2=field_x(1:401,1:181,1:181,1);% 860 nm data
ez2=field_x(1:401,1:181,1:181,1);% 860 nm data
%wavelength_2 = 3e8/field_f(1)
e22=abs(ex2).^2+abs(ey2).^2+abs(ez2).^2;

sim_origin = [170, 0, 0]*10^-9;% for nanosphere dimers different origin!!
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );
%dist1=sqrt((xmg+offset).^2+ymg.^2+zmg.^2);
%dist2=sqrt((xmg-offset).^2+ymg.^2+zmg.^2);
dist1=sqrt((xmg-offset).^2+ymg.^2+zmg.^2);
dist2=sqrt((xmg+offset).^2+ymg.^2+zmg.^2);

e2(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;
e21(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;
e22(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;

e2 = e21.*e22;

sum(sum(sum(e2)))

%% Smooth nanoshell dimer (separate E at 785 nm and 860 nm)
clear all
label1 = 'Smooth NS_R92_R169_dimer_sep_1nm_Air_3D_785_';
label2 = 'Smooth NS_R92_R169_dimer_sep_1nm_Air_3D_860_';

'Dont forget to change offset!!!'

offset = ((338+1)*10^-9)/2;
dist_min = 92e-9; %define E2 calculation region 
dist_max = 171e-9; % 170 for 1 nm over surface, 171 for 2 nm over surface 

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

sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;% for smooth nanoshell dimer R =169
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );
dist1=sqrt((xmg+offset).^2+ymg.^2+zmg.^2);
dist2=sqrt((xmg-offset).^2+ymg.^2+zmg.^2);

e2(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;
e21(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;
e22(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;

e2 = e21.*e22;
sum(sum(sum(e2)))

%% Small Au bead (r=13) single*2 
clear all
label1 = 'AuSphere_R13_dimer_sep50nm_Air_3D_785_860';

'Dont forget to change offset!!!'

offset = ((26+50)*10^-9)/2;
dist_min = 13e-9; %define E2 calculation region 
dist_max = 15e-9; % 14 for 1 nm over surface, 15 for 2 nm over surface 

load(sprintf('%sfieldx.mat', label1));
load(sprintf('%sfieldy.mat', label1));
load(sprintf('%sfieldz.mat', label1));
load(sprintf('%spos.mat', label1));

ex1=field_x(1:121,1:81,1:81,2);% 785 nm data
ey1=field_y(1:121,1:81,1:81,2);% 785 nm data
ez1=field_z(1:121,1:81,1:81,2);% 785 nm data
%wavelength_1 = 3e8/field_f(2)
e21=abs(ex1).^2+abs(ey1).^2+abs(ez1).^2;

ex2=field_x(1:121,1:81,1:81,1);% 860 nm data
ey2=field_x(1:121,1:81,1:81,1);% 860 nm data
ez2=field_x(1:121,1:81,1:81,1);% 860 nm data
%wavelength_2 = 3e8/field_f(1)
e22=abs(ex2).^2+abs(ey2).^2+abs(ez2).^2;

sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );
dist1=sqrt((xmg+offset).^2+ymg.^2+zmg.^2);
dist2=sqrt((xmg-offset).^2+ymg.^2+zmg.^2);

e2(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;
e21(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;
e22(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))))=0;

e2 = e21.*e22;

sum(sum(sum(e2)))

%% Small Au bead (r=13) trimer*2 
clear all
label1 = 'AuSphere_R13_trimer_sep50nm_Air_3D_785_860';

'Dont forget to change offset!!!'

offset_1 = ((26+50)*10^-9)/2;
offset_2 = ((28*2+26+50)*10^-9)/2;
offset_3 = (offset_1 + offset_2)/2;
dist_min = 13e-9; %define E2 calculation region 
dist_max = 15e-9; % 14 for 1 nm over surface, 15 for 2 nm over surface 

load(sprintf('%sfieldx.mat', label1));
load(sprintf('%sfieldy.mat', label1));
load(sprintf('%sfieldz.mat', label1));
load(sprintf('%spos.mat', label1));

ex1=field_x(1:221,1:101,1:81,2);% 785 nm data
ey1=field_y(1:221,1:101,1:81,2);% 785 nm data
ez1=field_z(1:221,1:101,1:81,2);% 785 nm data
%wavelength_1 = 3e8/field_f(2)
e21=abs(ex1).^2+abs(ey1).^2+abs(ez1).^2;

ex2=field_x(1:221,1:101,1:81,1);% 860 nm data
ey2=field_x(1:221,1:101,1:81,1);% 860 nm data
ez2=field_x(1:221,1:101,1:81,1);% 860 nm data
%wavelength_2 = 3e8/field_f(1)
e22=abs(ex2).^2+abs(ey2).^2+abs(ez2).^2;

sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );
dist1=sqrt((xmg + offset_1).^2 + (ymg + 8.08*10^-9).^2 + zmg.^2);
dist2=sqrt((xmg - offset_1).^2 + (ymg + 8.08*10^-9).^2 + zmg.^2);
dist3=sqrt((xmg + offset_2).^2 + (ymg + 8.08*10^-9).^2 + zmg.^2);
dist4=sqrt((xmg - offset_2).^2 + (ymg + 8.08*10^-9).^2 + zmg.^2);
dist5=sqrt((xmg + offset_3).^2 + (ymg - 16.16*10^-9).^2 + zmg.^2);
dist6=sqrt((xmg - offset_3).^2 + (ymg - 16.16*10^-9).^2 + zmg.^2);

e21(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))|((dist3<=dist_max)&(dist3>dist_min))|((dist4<=dist_max)&(dist4>dist_min))|((dist5<=dist_max)&(dist5>dist_min))|((dist6<=dist_max)&(dist6>dist_min))))=0;
e22(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))|((dist3<=dist_max)&(dist3>dist_min))|((dist4<=dist_max)&(dist4>dist_min))|((dist5<=dist_max)&(dist5>dist_min))|((dist6<=dist_max)&(dist6>dist_min))))=0;
e2(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))|((dist3<=dist_max)&(dist3>dist_min))|((dist4<=dist_max)&(dist4>dist_min))|((dist5<=dist_max)&(dist5>dist_min))|((dist6<=dist_max)&(dist6>dist_min))))=0;

e2 = e21.*e22;

sum(sum(sum(e2)))

%% Small Au bead (r=13) 4sphere*2 
clear all
label1 = 'AuSphere_R13_4sphere_sep50nm_Air_3D_785_860';

'Dont forget to change offset!!!'

offset_1 = ((26+50)*10^-9)/2;
offset_2 = ((28*2+26+50)*10^-9)/2;
dist_min = 13e-9; %define E2 calculation region 
dist_max = 15e-9; % 14 for 1 nm over surface, 15 for 2 nm over surface  

load(sprintf('%sfieldx.mat', label1));
load(sprintf('%sfieldy.mat', label1));
load(sprintf('%sfieldz.mat', label1));
load(sprintf('%spos.mat', label1));

ex1=field_x(1:221,1:101,1:81,2);% 785 nm data
ey1=field_y(1:221,1:101,1:81,2);% 785 nm data
ez1=field_z(1:221,1:101,1:81,2);% 785 nm data
%wavelength_1 = 3e8/field_f(2)
e21=abs(ex1).^2+abs(ey1).^2+abs(ez1).^2;

ex2=field_x(1:221,1:101,1:81,1);% 860 nm data
ey2=field_y(1:221,1:101,1:81,1);% 860 nm data
ez2=field_z(1:221,1:101,1:81,1);% 860 nm data
%wavelength_2 = 3e8/field_f(1)
e22=abs(ex2).^2+abs(ey2).^2+abs(ez2).^2;

sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );
dist1=sqrt((xmg + offset_1).^2 + (ymg + 14*10^-9).^2 + zmg.^2);
dist2=sqrt((xmg - offset_1).^2 + (ymg + 14*10^-9).^2 + zmg.^2);
dist3=sqrt((xmg + offset_2).^2 + (ymg + 14*10^-9).^2 + zmg.^2);
dist4=sqrt((xmg - offset_2).^2 + (ymg + 14*10^-9).^2 + zmg.^2);
dist5=sqrt((xmg + offset_1).^2 + (ymg - 14*10^-9).^2 + zmg.^2);
dist6=sqrt((xmg - offset_1).^2 + (ymg - 14*10^-9).^2 + zmg.^2);
dist7=sqrt((xmg + offset_2).^2 + (ymg - 14*10^-9).^2 + zmg.^2);
dist8=sqrt((xmg - offset_2).^2 + (ymg - 14*10^-9).^2 + zmg.^2);

e21(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))|((dist3<=dist_max)&(dist3>dist_min))|((dist4<=dist_max)&(dist4>dist_min))|((dist5<=dist_max)&(dist5>dist_min))|((dist6<=dist_max)&(dist6>dist_min))|((dist7<=dist_max)&(dist7>dist_min))|((dist8<=dist_max)&(dist8>dist_min))))=0;
e22(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))|((dist3<=dist_max)&(dist3>dist_min))|((dist4<=dist_max)&(dist4>dist_min))|((dist5<=dist_max)&(dist5>dist_min))|((dist6<=dist_max)&(dist6>dist_min))|((dist7<=dist_max)&(dist7>dist_min))|((dist8<=dist_max)&(dist8>dist_min))))=0;
e2(~(((dist1<=dist_max)&(dist1>dist_min))|((dist2<=dist_max)&(dist2>dist_min))|((dist3<=dist_max)&(dist3>dist_min))|((dist4<=dist_max)&(dist4>dist_min))|((dist5<=dist_max)&(dist5>dist_min))|((dist6<=dist_max)&(dist6>dist_min))|((dist7<=dist_max)&(dist7>dist_min))|((dist8<=dist_max)&(dist8>dist_min))))=0;

e2 = e21.*e22;

sum(sum(sum(e2)))

%% Spherical Au NP 
clear all
label1 = 'AuSphere_R169_Air_3D_785_860';

'Dont forget to change offset!!!'

dist_min = 92e-9; %define E2 calculation region 
dist_max = 170e-9; % 170 for 1 nm over surface, 171 for 2 nm over surface 

load(sprintf('%sfieldx.mat', label1));
load(sprintf('%sfieldy.mat', label1));
load(sprintf('%sfieldz.mat', label1));
load(sprintf('%spos.mat', label1));

ex1=field_x(1:181,1:181,1:181,2);% 785 nm data
ey1=field_y(1:181,1:181,1:181,2);% 785 nm data
ez1=field_z(1:181,1:181,1:181,2);% 785 nm data
%wavelength_1 = 3e8/field_f(2)
e21=abs(ex1).^2+abs(ey1).^2+abs(ez1).^2;

ex2=field_x(1:181,1:181,1:181,1);% 860 nm data
ey2=field_x(1:181,1:181,1:181,1);% 860 nm data
ez2=field_x(1:181,1:181,1:181,1);% 860 nm data
%wavelength_2 = 3e8/field_f(1)
e22=abs(ex2).^2+abs(ey2).^2+abs(ez2).^2;

sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;% for nanosphere dimers
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );
dist=sqrt(xmg.^2+ymg.^2+zmg.^2);

e21(~((dist<=dist_max)&(dist>dist_min)))=0;
e22(~((dist<=dist_max)&(dist>dist_min)))=0;
e2(~((dist<=dist_max)&(dist>dist_min)))=0;

e2 = e21.*e22;

sum(sum(sum(e2)))