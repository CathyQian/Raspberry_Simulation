clear all
label = 'BDAC_N800_R92_r14_ps1_dimer_sep_1nm_Air_3D_785_';
savemap = 'TBDAC_N800_R92_r14_ps1_dimer_sep_1nm_Air_3D_785_E_particle';

cutoff = 2;
circlesize = 16;
fontsize = 45;

load(sprintf('%sfieldx.mat', label));
load(sprintf('%sfieldy.mat', label));
load(sprintf('%sfieldz.mat', label));
load(sprintf('%spos.mat', label));
%load(sprintf('%sindex.mat', label));

ex=field_x;
ey=field_y;
ez=field_z;
e2=abs(ex).^2+abs(ey).^2+abs(ez).^2;
sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;% for 184 nm core and
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );

xp=reshape(xmg,size(xmg,1)*size(xmg,2)*size(xmg,3),1,1);
yp=reshape(ymg,size(ymg,1)*size(ymg,2)*size(ymg,3),1,1);
zp=reshape(zmg,size(zmg,1)*size(zmg,2)*size(zmg,3),1,1);
%real_index=reshape(real(ind_ix),size(xmg,1)*size(xmg,2)*size(xmg,3),1,1);

x1=xp;
y1=yp;
z1=zp;
%
e2lin = reshape(e2,size(e2,1)*size(e2,2)*size(e2,3),1,1);
%h1=figure;
e2linlog = log10(e2lin);
%hist(e2linlog); %Use this line to view histogram for the log scale E field
%intensity. Use it to define cutoff
%
e2linlogsel = e2linlog;

dist = sqrt(xp.^2+yp.^2+zp.^2);
%xp(dist<48e-9)=[];
%yp(dist<48e-9)=[];
%zp(dist<48e-9)=[];
%real_index(dist<48e-9)=[];
%e2linlogsel(dist<48e-9)=[];

xp(zp<0e-8)=[];
yp(zp<0e-8)=[];
%real_index(zp<0e-8)=[];
e2linlogsel(zp<0e-8)=[];
zp(zp<0e-8)=[];

% x(y<0)=[];
% z(y<0)=[];
% real_index(y<0)=[];
% e2linlogsel(y<0)=[];
% y(y<0)=[];
% 
% 
% y(x<0)=[];
% z(x<0)=[];
% real_index(x<0)=[];
% e2linlogsel(x<0)=[];
% x(x<0)=[];

%xp(real_index<1.4)=[];
%yp(real_index<1.4)=[];
%zp(real_index<1.4)=[];
%e2linlogsel(real_index<1.4)=[];

xp(e2linlogsel<cutoff)=[];
yp(e2linlogsel<cutoff)=[];
zp(e2linlogsel<cutoff)=[];

e2linlogsel(e2linlogsel<cutoff)=[];

NumOfDataPts = length(e2linlogsel)

%view(45,45)
%%
%colorvals = (e2linlogsel-min(e2linlogsel));
%colorvals = colorvals/max(colorvals);

colormap(jet)
cm=colormap;
scale = linspace(cutoff,5.2,64);

for i=1:length(e2linlogsel)
cmp=e2linlogsel(i)>=scale;
cnt=sum(cmp);
cmarker(i,:)=cm(cnt,:);
end
%h2=scatter3(x,y,z,circlesize,e2linlogsel,'fill');

%hp = figure;
hold on
for i = 1:length(xp)
    plot3(xp(i)*1e9,yp(i)*1e9,zp(i)*1e9, '.', 'markersize', 10, 'MarkerEdgeColor',cmarker(i,:));
end

xlabel('X (nm)', 'fontSize', fontsize, 'fontWeight', 'bold');
ylabel('Y (nm)', 'fontSize', fontsize, 'fontWeight','bold');
set(gca,'fontsize',fontsize, 'fontWeight', 'bold','LineWidth',2);

%saveas(hp,savemap,'tif'); 
%saveas(hp,savemap,'fig'); 
grid off
%view([45,45])
%saveas(h2,sprintf('%s_Esq_3d_plot.jpg',label));

%%
%load('dt_0.001000 damp_30.000000 N_800 k=80.000000 br=14.000000 cr=92.000000 var=1.mat');
load('Two_BDAC_offset_dist_330_184core_GNC_sep_1nm_copy');
[Xs,Ys,Zs]=sphere(20);
%hf=figure;
%ha=axes;
hold on
for i=1:1600
xs1=(Xs*16e-9)+x(i,1)*1e-9;
ys1=(Ys*16e-9)+x(i,2)*1e-9;
zs1=(Zs*16e-9)+x(i,3)*1e-9;

    for j=1:21
        for k=1:21
            distZ(j,k)= abs(zs1(j,k)-x(i,3)*1e-9);
        end
    end
    h(i)=surf(xs1,ys1,zs1,distZ);
end


for i=1:1600
set(h(i),'edgecolor','none')
set(h(i),'facealpha',0.5)
end

colormap(summer);
