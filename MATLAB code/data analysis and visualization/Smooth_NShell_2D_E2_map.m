% generate 2D E2 at 785 nm or 860 nm map for smooth nanoshells
% try not to use cutoff in order to show all data points

clear all
label = 'Smooth NS_R92_R169_dimer_sep_5nm_Air_2D_E_';

cutoff = 0; % don't forget to change cutoff

load(sprintf('%sfieldx.mat', label));
load(sprintf('%sfieldy.mat', label));
load(sprintf('%sfieldz.mat', label));
load(sprintf('%spos.mat', label));

ex_785=field_x(: , : , : ,105);% 400~1200 nm range, 201 data points, 785 nm is 105th and 860 nm is 87th data point
ey_785=field_y(: , : , : ,105);
ez_785=field_z(: , : , : ,105);
e2_785=abs(ex_785).^2+abs(ey_785).^2+abs(ez_785).^2;

ex_860=field_x(: , : , : ,87);
ey_860=field_y(: , : , : ,87);
ez_860=field_z(: , : , : ,87);
e2_860=abs(ex_860).^2+abs(ey_860).^2+abs(ez_860).^2;

sim_origin = [-4.37751, 35.3597, -36.9736]*10^-9;% for smooth nanoshells R92_R169
field_xp = field_xp - sim_origin(1);
field_yp = field_yp - sim_origin(2);
field_zp = field_zp - sim_origin(3);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );% transfer vectors into arrays

xp=reshape(xmg,size(xmg,1)*size(xmg,2)*size(xmg,3),1,1);% change position matrix into 1D array
yp=reshape(ymg,size(ymg,1)*size(ymg,2)*size(ymg,3),1,1);
zp=reshape(zmg,size(zmg,1)*size(zmg,2)*size(zmg,3),1,1);

e2_785lin = reshape(e2_785,size(e2_785,1)*size(e2_785,2)*size(e2_785,3),1,1);
e2_860lin = reshape(e2_860,size(e2_860,1)*size(e2_860,2)*size(e2_860,3),1,1);

e2_785linlogsel = log10(e2_785lin);
e2_860linlogsel = log10(e2_860lin);

e2_785linlogsel(e2_785linlogsel<cutoff)=[];
e2_860linlogsel(e2_860linlogsel<cutoff)=[];
xp(e2_785linlogsel<cutoff)=[];
yp(e2_785linlogsel<cutoff)=[];
zp(e2_785linlogsel<cutoff)=[];


NumOfDataPts_785 = length(e2_785linlogsel)
NumOfDataPts_860 = length(e2_860linlogsel)

%%

colormap(jet) % set color coding
cm=colormap;
%scale = linspace(min(e2_785linlogsel),max(e2_785linlogsel),64);
scale = linspace(cutoff,max(e2_785linlogsel),64);

for i=1:length(e2_785linlogsel) % set color for each log10(E2) value
cmp=e2_785linlogsel(i)>=scale;
cnt=sum(cmp);
cmarker(i,:)=cm(cnt,:);
end

%circlesize = 10;
for i = 1:length(xp)
      plot3(xp(i),yp(i),zp(i), '.', 'markersize', 10, 'MarkerEdgeColor',cmarker(i,:));
end
grid off
%%

for i=1:length(e2_860linlogsel)
cmp=e2_860linlogsel(i)>=scale;
cnt=sum(cmp);
cmarker(i,:)=cm(cnt,:);
end


for i = 1:length(xp)
    scatter(xp(i),yp(i),circlesize,cmarker(i,:));
end
grid off

