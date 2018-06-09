clear all
label1 = 'BDAC_N800_R92_r14_ps1_dimer_sep_1nm__Air_3D_785_p_';
label2 = 'BDAC_N800_R92_r14_ps1_dimer_sep_1nm__Air_3D_860_p_';


load(sprintf('%sfieldx.mat', label1));
load(sprintf('%sfieldy.mat', label1));
load(sprintf('%sfieldz.mat', label1));

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

e4 = e21.*e22;

sum(sum(sum(e4)))