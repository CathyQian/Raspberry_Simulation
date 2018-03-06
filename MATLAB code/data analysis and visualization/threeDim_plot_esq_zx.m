%% calcalute the E2 at 785 nm for individual GNCs or GNC dimers

 clear all
label = 'BDAC_N800_R92_r14_ps1_Air_3D_785_';
savemap = 'BDAC_N800_R92_r14_ps1_Air_3D_785_';
%savehis = 'BDAC_N800_R92_r14_ps1_dimer_sep_-2nm__Air_860_E2_his';

fontsize = 45;
cutoff = 2; % don't forget to change cutoff
%circlesize = 16;

load(sprintf('%sfieldx.mat', label));
load(sprintf('%sfieldy.mat', label));
load(sprintf('%sfieldz.mat', label));
load(sprintf('%spos.mat', label));
%load(sprintf('%sindex.mat', label));

ex=field_x;
ey=field_y;
ez=field_z;
e2=abs(ex).^2+abs(ey).^2+abs(ez).^2;
%sim_origin = [20.6017, 32.8531, 45.8893]*10^-9;% for 96 nm core
sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;% for 184 nm core and
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );

xp=reshape(xmg,size(xmg,1)*size(xmg,2)*size(xmg,3),1,1);
yp=reshape(ymg,size(ymg,1)*size(ymg,2)*size(ymg,3),1,1);
zp=reshape(zmg,size(zmg,1)*size(zmg,2)*size(zmg,3),1,1);
%real_index=reshape(real(ind_ix),size(xmg,1)*size(xmg,2)*size(xmg,3),1,1);
%no index data

x1=xp;
y1=yp;
z1=zp;
%
e2lin = reshape(e2,size(e2,1)*size(e2,2)*size(e2,3),1,1);
%h1=figure;% figure handle,used for /processing saving figure
e2linlog = log10(e2lin);
%hist(e2linlog); %Use this line to view histogram for the log scale E field
%intensity. Use it to define cutoff

%%xlim([-4,6]);
%xlabel('Log_{10}(E^{2} at 860 nm)', 'fontSize', fontsize, 'fontWeight', 'bold');
%ylabel('Counts', 'fontSize', fontsize, 'fontWeight','bold');
%set(gca,'fontsize',fontsize, 'fontWeight', 'bold','LineWidth',2);
%h2 = findobj(gca,'Type','patch');
%set(h2,'FaceColor','k','EdgeColor','w');

%saveas(h1, savehis, 'tif');
%saveas(h1, savehis, 'fig');

%%
e2linlogsel = e2linlog;

dist = sqrt(xp.^2+yp.^2+zp.^2);
%xp(dist<48e-9)=[];% inside of PS core, there's no field
%yp(dist<48e-9)=[];
%zp(dist<48e-9)=[];
%real_index(dist<48e-9)=[];
%e2linlogsel(dist<48e-9)=[];%inside of PS core, the E2=0

xp(zp<0e-8)=[];% only draw half of the image (zp>0 region)
yp(zp<0e-8)=[];
%real_index(zp<0e-8)=[];
e2linlogsel(zp<0e-8)=[];
zp(zp<0e-8)=[];

xp(e2linlogsel<cutoff)=[];
yp(e2linlogsel<cutoff)=[];
zp(e2linlogsel<cutoff)=[];

e2linlogsel(e2linlogsel<cutoff)=[];

NumOfDataPts = length(e2linlogsel)

%%
%colorvals = (e2linlogsel-min(e2linlogsel));
%colorvals = colorvals/max(colorvals);

colormap(jet) %
cm=colormap;
scale = linspace(2,5.2,64);% don't forget to change, colorcoding range 2 to 5.2


for i=1:length(e2linlogsel)
cmp=e2linlogsel(i)>=scale;
cnt=sum(cmp);
cmarker(i,:)=cm(cnt,:);
end

%h2=scatter3(xp,yp,zp,circlesize,e2linlogsel,'fill');
%hold all
hp = figure;
%axis equal;
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
%This script creates 2 identical plot for colorbars: 1 with ticks and
%numbers as reference, the other without for actual use. Need to change the
%minimum and maximum of the scale, and number of ticks wanted
clear all
 
min = cutoff;
max = 6.2;% for 94.5 nm core, need to modify for 184 nm core
number_of_ticks = 5;

scale = linspace(min,max,number_of_ticks);
mscale = repmat([1:64],8,1);

hcf2 = figure; % with scale bar
hca2 = axes;
hcb2 = pcolor(mscale');
colormap(jet);
set(hcb2,'facecolor','interp');
set(hca2,'xtick',[]);
set(hca2,'ytick',linspace(1,64,number_of_ticks));
set(hca2,'yticklabel',scale);
set(hcb2,'edgecolor','none');
set(hca2,'yaxislocation','r');
set(hca2,'fontsize',fontsize);
set(hca2,'linewidth',0.01);
axis image
saveas(hcb2, 'colorbar_core94p5_dis_-1nm', 'tif');
saveas(hcb2, 'colorbar_core94p5_dis_-1nm', 'fig');
%%
