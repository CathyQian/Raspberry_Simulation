% calcalute the E2 at 785 nm for individual GNCs or GNC dimers

clear all
label = 'Smooth NS_R92_R169_dimer_sep_1nm_Air_3D_860_';

%cutoff = -2; % don't forget to change cutoff

load(sprintf('%sfieldx.mat', label));
load(sprintf('%sfieldy.mat', label));
load(sprintf('%sfieldz.mat', label));
load(sprintf('%spos.mat', label));

ex=field_x;
ey=field_y;
ez=field_z;
e2=abs(ex).^2+abs(ey).^2+abs(ez).^2;
sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;% for smooth nanoshells R92_R169
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
%h1=figure;
e2linlog = log10(e2lin);
%hist(e2linlog); %Use this line to view histogram for the log scale E field
%intensity. Use it to define cutoff
e2linlogsel = e2linlog;

%dist = sqrt(xp.^2+yp.^2+zp.^2);
%xp(dist<92e-9)=[];% inside of PS core, there's no field
%yp(dist<92e-9)=[];
%zp(dist<92e-9)=[]; doesn't work for dimers
%real_index(dist<48e-9)=[];
%e2linlogsel(dist<92e-9)=[];%inside of PS core, the E2=0

xp(zp>0.5e-9)=[]; %only plot data on -2<z<2 nm range since mesh = 2nm
xp(zp<-0.5e-9)=[];
yp(zp>0.5e-9)=[];
yp(zp<-0.5e-9)=[];
e2linlogsel(zp>0.5e-9)=[];
e2linlogsel(zp<-0.5e-9)=[];


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

%xp(real_index<1.4)=[];% only draw area with PS.
%yp(real_index<1.4)=[];% E2 inside of Au is zero.
%zp(real_index<1.4)=[];
%e2linlogsel(real_index<1.4)=[];

%xp(e2linlogsel<cutoff)=[];
%yp(e2linlogsel<cutoff)=[];
%zp(e2linlogsel<cutoff)=[];

%e2linlogsel(e2linlogsel<cutoff)=[];

NumOfDataPts = length(e2linlogsel)

%%
%colorvals = (e2linlogsel-min(e2linlogsel));
%colorvals = colorvals/max(colorvals);

colormap(jet) %
cm=colormap;
scale = linspace(min(e2linlogsel),max(e2linlogsel),64);

for i=1:length(e2linlogsel)
cmp=e2linlogsel(i)>=scale;
cnt=sum(cmp);
cmarker(i,:)=cm(cnt,:);
end

circlesize = 10;
%h2=scatter3(xp,yp,zp,circlesize,e2linlogsel,'fill');
hold on
for i = 1:length(xp)
     plot3(xp(i),yp(i),zp(i), '.', 'markersize', 10, 'MarkerEdgeColor',cmarker(i,:));
end
grid off
%view([45,45])
%saveas(h2,sprintf('%s_Esq_3d_plot.jpg',label));
