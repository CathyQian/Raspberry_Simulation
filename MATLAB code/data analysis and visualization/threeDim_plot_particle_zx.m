%clear all
load('dt_0.001000 damp_30.000000 N_800 k=80.000000 br=14.000000 cr=92.000000 var=1.mat');
[Xs,Ys,Zs]=sphere(20);
%hf=figure;
%ha=axes;
hold on
for i=1:800
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


for i=1:800
set(h(i),'edgecolor','none')
set(h(i),'facealpha',0.5)
end

colormap(summer);