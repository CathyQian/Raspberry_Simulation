clear all
 
min = 2;
max = 6.2;% for 94.5 nm core, need to modify for 184 nm core
number_of_ticks = 5;
fontsize = 45;

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
