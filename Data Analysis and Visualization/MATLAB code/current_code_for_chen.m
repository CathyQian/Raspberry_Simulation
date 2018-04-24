% Makes the SERS paper figures



%% Calculate dipole and quadrapole modes
clear all;
close all;
eps0 = 8.85418782e-12;

%label = 'tj_final_bc_air_';
label = 'Mode_ana_test_single_gold_br25_air_ms1_';


xfile = sprintf('%sfieldx.mat', label);
yfile = sprintf('%sfieldy.mat', label);
zfile = sprintf('%sfieldz.mat', label);

ffile = sprintf('%sfreq.mat', label);
pfile = sprintf('%spos.mat', label);       

ifile = sprintf('%sindex.mat', label);

load(xfile);
load(yfile);
load(zfile);

load(ffile);
load(pfile);
load(ifile);


%field_xp = field_xp -( -6.40585e-9);
%field_yp = field_yp -( 40.1131e-9);
%field_zp = field_zp -( -3.68297e-9);
field_xp = field_xp -(0);
field_yp = field_yp -(0);
field_zp = field_zp -(0);
[xmg, ymg, zmg] = ndgrid( field_xp, field_yp, field_zp );

% get position vectors
f = field_f;

% get the complex index and fields. Assume it's isotropic
epsilon_x = eps0*ind_ix.^2;
epsilon_y = eps0*ind_ix.^2;
epsilon_z = eps0*ind_ix.^2;

Ex = field_x;
Ey = field_y;
Ez = field_z;

Dx = epsilon_x .* Ex;
Dy = epsilon_y .* Ey;
Dz = epsilon_z .* Ez;
%D0x = eps0 * Ex;
%D0y = eps0 * Ey;
%D0z = eps0 * Ez;
%
Jx = zeros(size(Dx));
Jy = zeros(size(Dx));
Jz = zeros(size(Dx));
% calculate the current(bound)
for k = 1:length(f)
    Jx(:,:,:,k) = 1i*2*pi*f(k) * (eps0*Ex(:,:,:,k) - Dx(:,:,:,k));
    Jy(:,:,:,k) = 1i*2*pi*f(k) * (eps0*Ey(:,:,:,k) - Dy(:,:,:,k));
    Jz(:,:,:,k) = 1i*2*pi*f(k) * (eps0*Ez(:,:,:,k) - Dz(:,:,:,k));
end
%
field_int = ( abs(Ex( :, :, :,:)).^2 ...
+abs(Ey( :, :, :,:)).^2 ...
+abs(Ez( :, :, :,:)).^2 );


%%


% section for writing out color maps with matlab
close all;

[s1, s2, s3, s4] = size(field_int);
vol = s3;
mx = 50;

% Set of wavelengths to use for spikies
% d = [30, 784; 69,634; 95, 562; 85, 588; 89, 576; 50,700; 57, 672; 61, 660;48, 707 ];
% set to use for QA Spheres:
d = [ 85,588];

for i = 1:length(d)
    k=d(i,1);
    filename = sprintf('%s_E2int_%dmax_%d.tiff', label,mx, d(i,2) );
    
    figureHandle = figure; imagesc( [118:-2:-118], [-118:2:118],sum(field_int( :, :, :, k )./vol,3)');
    colorbar;
    
    
    caxis( [0, mx] );
    CM=get(gcf,'colormap');
    colormap(flipud(CM))
            

    xlhand = get(gca,'xlabel');
    set(xlhand,'string','x(nm)','fontsize',22)
    xlhand = get(gca,'ylabel');
    set(xlhand,'string', 'y(nm)', 'fontsize',22)
    set(gca,'fontsize',22)
    set(gcf,'PaperPositionMode','auto')
    
    %axis image
    axis xy
    %set(figureHandle, 'Position', [0 0 900 900])
    
    set( gca, 'XTick', [-100, -50, 0, 50, 100] );
    set( gca, 'YTick', [-100, -50, 0, 50, 100] );
    set(figureHandle,'color','w');
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    f=getframe(figureHandle);
    
    imwrite(f.cdata, filename, 'Resolution', 300);
end

%% Code for imaging the currents(outputs a text file)
[s1, s2, s3, s4] = size(field_int);
vol = s3;
mx = 50;

% Set of wavelengths to use for spikies
% d = [30, 784; 69,634; 95, 562; 85, 588; 89, 576; 50,700; 57, 672; 61, 660;48, 707 ];
% set to use for QA Spheres:
d = [ 85,588];
for i = 1:length(d)
    k=d(i,1);
    filename = sprintf('%s_J_%d.txt', label, d(i,2) );    

    Jxs = real(Jx(:,:,:,k));
    Jys = real(Jy(:,:,:,k));    

    %project onto XY plane
    Jxs_proj = squeeze(sum( Jxs, 3));
    Jys_proj = squeeze(sum( Jys, 3));
    Jmag = sqrt( Jxs_proj.^2 + Jys_proj.^2 );
    
    lin_ind_p = 1:(length(Jmag)^2);
    save_id = rand(length(Jmag).^2,1)<0.2;
    lin_ind = lin_ind_p(save_id);
    [xvals,yvals,zvals]= ind2sub(size(Jmag),lin_ind);
    Jnorm = 1./Jmag(lin_ind);
    Jnorm(isnan(Jnorm)) = 0;
    Jxout = (Jxs_proj(lin_ind).*Jnorm)';
    Jyout = (Jys_proj(lin_ind).*Jnorm)';
    Jxout(isnan(Jxout)) = 0;
    Jyout(isnan(Jyout)) = 0;
    out = [ ((-xvals*2)+118)', ((yvals*2)-118)', atan2(Jyout, -Jxout), (abs( Jxout) + abs(Jyout) > 1e-5)*12 ];
    save( filename, '-ascii', 'out' );
end



%%
clear P
% Find the dipole
for k = 1:length(f)
    P(k, 1:3) = 1*10^-27*1i./f(k)/2/3.14159*[ sum(sum(sum(Jx(:,:,:,k)))), sum(sum(sum(Jy(:,:,:,k)))),sum(sum(sum(Jz(:,:,:,k))))];
end
wl = 3e8./f;

Pabs = abs(P);

Pnorm = sqrt( Pabs(:,1).^2+Pabs(:,2).^2+Pabs(:,3).^2);

figure;
plot( wl, Pnorm )

%
figure;
hold on;
plot( wl, Pabs(:,1),'x');
plot( wl, Pabs(:,2),'o');
plot( wl, Pabs(:,3),'s');

%
clear Q
% Find the quadrapole
for k = 1:length(f)
    Q(k, 1, 1) = 1i./f(k)/2/3.14159*[ sum(sum(sum(Jx(:,:,:,k).*xmg))) + sum(sum(sum(Jx(:,:,:,k).*xmg)))]*1*10^-27;
    Q(k, 1, 2) = 1i./f(k)/2/3.14159*[ sum(sum(sum(Jx(:,:,:,k).*ymg))) + sum(sum(sum(Jy(:,:,:,k).*xmg)))]*1*10^-27;
    Q(k, 1, 3) = 1i./f(k)/2/3.14159*[ sum(sum(sum(Jx(:,:,:,k).*zmg))) + sum(sum(sum(Jz(:,:,:,k).*xmg)))]*1*10^-27;    

    Q(k, 2, 1) = 1i./f(k)/2/3.14159*[ sum(sum(sum(Jy(:,:,:,k).*xmg))) + sum(sum(sum(Jx(:,:,:,k).*ymg)))]*1*10^-27;
    Q(k, 2, 2) = 1i./f(k)/2/3.14159*[ sum(sum(sum(Jy(:,:,:,k).*ymg))) + sum(sum(sum(Jy(:,:,:,k).*ymg)))]*1*10^-27;
    Q(k, 2, 3) = 1i./f(k)/2/3.14159*[ sum(sum(sum(Jy(:,:,:,k).*zmg))) + sum(sum(sum(Jz(:,:,:,k).*ymg)))]*1*10^-27;        
    
    Q(k, 3, 1) = 1i./f(k)/2/3.14159*[ sum(sum(sum(Jz(:,:,:,k).*xmg))) + sum(sum(sum(Jx(:,:,:,k).*zmg)))]*1*10^-27;
    Q(k, 3, 2) = 1i./f(k)/2/3.14159*[ sum(sum(sum(Jz(:,:,:,k).*ymg))) + sum(sum(sum(Jy(:,:,:,k).*zmg)))]*1*10^-27;
    Q(k, 3, 3) = 1i./f(k)/2/3.14159*[ sum(sum(sum(Jz(:,:,:,k).*zmg))) + sum(sum(sum(Jz(:,:,:,k).*zmg)))]*1*10^-27;        
    
end
%%
close all;
Qabs = abs(Q);
for( k = 1:length(f) )
    Qnorm(k) = norm(squeeze(Qabs(k,:,:)));
end


figure;
plot( wl, Qnorm )

figure; hold on;
plot( wl, Qabs(:,1,1),'sr');
plot( wl, Qabs(:,1,2),'-xk');
plot( wl, Qabs(:,1,3),'--ob');
%plot( wl, Qabs(:,2,1),'o--r');
%plot( wl, Qabs(:,2,2),'o--k');
%plot( wl, Qabs(:,2,3),'o-.b');
%plot( wl, Qabs(:,3,1),'-.r');
%plot( wl, Qabs(:,3,2),'-.k');
%plot( wl, Qabs(:,3,3),'-.b');
legend('xx', 'xy', 'xz' );
%%
%
clear m
% Find the magnetic dipole
for k = 1:length(f)
    m(k, 1:3) = 1/2*[ sum(sum(sum(Jz(:,:,:,k).*ymg - Jy(:,:,:,k).*zmg)))*1.0*10^-27,...
        sum(sum(sum(Jx(:,:,:,k).*zmg - Jz(:,:,:,k).*xmg)))*1.0*10^-27, ...
        sum(sum(sum(Jy(:,:,:,k).*xmg - Jx(:,:,:,k).*ymg)))*1.0*10^-27];
end
%%
figure;
hold on;
plot( wl, abs(m(:,1)),'o--k');
plot( wl, abs(m(:,2)),'x--r');
plot( wl, abs(m(:,3)),'-b');
legend('x', 'y', 'z');
%%
output_data = [wl, Qabs(:,1,2) ];
fileout = sprintf('%sQXYResult.txt', label);
save(fileout, '-ascii', 'output_data');

output_data = [wl, Qabs(:,2,2) ];
fileout = sprintf('%sQYYResult.txt', label);
save(fileout, '-ascii', 'output_data');

output_data = [wl, Pabs ];
fileout = sprintf('%sPResult.txt', label);
save(fileout, '-ascii', 'output_data');



%% Calculate E^2
clear all;

%label = 'tj_final_bc_air_';
label = 'r4_final_bc_air_';


xfile = sprintf('%sfieldx.mat', label);
yfile = sprintf('%sfieldy.mat', label);
zfile = sprintf('%sfieldz.mat', label);

ffile = sprintf('%sfreq.mat', label);
pfile = sprintf('%spos.mat', label);         % Since the monitor doesn't move, these are all the same
ifile = sprintf('%sindex.mat', label);

load(xfile);
load(yfile);
load(zfile);

load(ffile);
load(pfile);
load(ifile);

full_field = cat(5, field_x, field_y, field_z);
size(full_field)

% Combine these pieces of data

field_magnitude = sqrt( abs(full_field( :, :, :,:,1)).^2 ...
    +abs(full_field( :, :, :,:,2)).^2 ...
    +abs(full_field( :, :, :,:,3)).^2 );


wl2 = 3e8./field_f.*1e9;

result = squeeze(sum(sum(sum(field_magnitude))));
size(squeeze(result))

field_xp = field_xp -( -6.40585e-9);
field_yp = field_yp -( 40.1131e-9);
field_zp = field_zp -( -3.68297e-9);
[fpx, fpy, fpz] = ndgrid( field_xp, field_yp, field_zp );

fcoords = [fpx(:), fpy(:), fpz(:)];
[A,B,C,D] = size( ind_ix );
if( D > 1 )
    ind_val = abs(ind_ix(:,:,:,floor(D/2)));
else
    ind_val = abs(ind_ix(:,:,:));
end
ind_x = ind_x -( -6.40585e-9 );
ind_y = ind_y -( 40.1131e-9);
ind_z = ind_z -( -3.68297e-9);


[x,y,z] = ndgrid( ind_x, ind_y, ind_z );
rad = sqrt( x.^2 + y.^2 + z.^2 );
overwrite_vals = rad < 54e-9;
ind_val(overwrite_vals) = 2;

iss = isosurface(x,y,z,abs(ind_val(:,:,:)), 1.4);

tic
DT = DelaunayTri(iss.vertices(:,1), iss.vertices(:,2), iss.vertices(:,3));
pid = nearestNeighbor(DT, fcoords);
toc

ncoords = iss.vertices(pid,:);

dist = sqrt( (ncoords(:,1) - fcoords(:,1)).^2 + (ncoords(:,2) - fcoords(:,2)).^2 + (ncoords(:,3) - fcoords(:,3)).^2 );
exterior = sqrt( (ncoords(:,1)).^2 + (ncoords(:,2)).^2 + (ncoords(:,3)).^2 ) < sqrt( (fcoords(:,1)).^2 + (fcoords(:,2)).^2 + (fcoords(:,3)).^2 ); % 1 or true if it's on the outside. This is fucking miraculous.
ext_scale = (exterior - 0.5)*2; % interior points are -1, exterior points are +1
int_scale = -1*ext_scale;
dist = int_scale.*dist;


cut_off1 = 0;
cut_off2 = -2.1;

for k = 1:1:length(field_f)
    fftemp = field_magnitude(:,:,:,k);
    summed_value(k) = sum( fftemp( dist < cut_off1*1e-9 & dist > cut_off2*1e-9));
    summed_value_E2_02(k) = sum( fftemp( dist < cut_off1*1e-9 & dist > cut_off2*1e-9).^2);
end
in_save = dist < cut_off1*1e-9 & dist > cut_off2*1e-9;
'Number Pixels:'
sum(in_save)

output_data = [wl2, summed_value_E2_02' ];
fileout = sprintf('%sE2Result.txt', label);
save(fileout, '-ascii', 'output_data');


j=1;    
for k = 1:1:length(wl2)
    wl_exc = wl2(k);
    wl_scat = 1./( 1./wl_exc - 1077*10^-7);
    if wl_scat > wl2(1)
        continue;
    end
    fftemp_exc = field_magnitude(:,:,:,find(wl2>wl_exc, 1, 'last'));
    fftemp_sca = field_magnitude(:,:,:,find(wl2>wl_scat, 1, 'last'));

    
    ef_1077(j) = sum( fftemp_exc( dist < cut_off1*1e-9 & dist > cut_off2*1e-9).^2 .* fftemp_sca( dist < cut_off1*1e-9 & dist > cut_off2*1e-9).^2);    
    ef_1077wl(j) = wl2(k);   
    j=j+1;
end
output_data = [ef_1077wl',  ef_1077' ];
fileout = sprintf('%sE4_1077.txt', label);
save(fileout, '-ascii', 'output_data');



j=1;    
for k = 1:1:length(wl2)
    wl_exc = wl2(k);
    wl_scat = 1./( 1./wl_exc - 1585*10^-7);
    if wl_scat > wl2(1)
        continue;
    end
 
    
    fftemp_exc = field_magnitude(:,:,:,find(wl2>wl_exc, 1, 'last'));
    fftemp_sca = field_magnitude(:,:,:,find(wl2>wl_scat, 1, 'last'));

    ef_1585(j) = sum( fftemp_exc( dist < cut_off1*1e-9 & dist > cut_off2*1e-9).^2 .* fftemp_sca( dist < cut_off1*1e-9 & dist > cut_off2*1e-9).^2);    
    ef_1585wl(j) = wl2(k);   
    j=j+1;
end

output_data = [ef_1585wl',  ef_1585' ];
fileout = sprintf('%sE4_1585.txt', label);
save(fileout, '-ascii', 'output_data');

%% This stuff just more or less gives us a record that stuff looked good
 
   

filename = sprintf('%s_3D_render.tiff', label );     
figureHandle = figure;
[A,B,C] = size(fpx);
dist3d = reshape(dist, [A,B,C] );
int_pts = reshape(ext_scale, [A,B,C] );
iss2 = isosurface( fpx,fpy,fpz, dist3d,-2e-9 );
p2 = patch(iss2);
set(p2, 'FaceColor', 'yellow', 'EdgeColor', 'none');
daspect([1 1 1])
view(3)
camlight; lighting phong
axis off;
f=getframe(figureHandle);
imwrite(f.cdata, filename, 'Resolution', 300);

filename = sprintf('%s_index_check60.tiff', label );     
figureHandle = figure;
imagesc( ind_val( :,  :,   60 )<1.4);
f=getframe(figureHandle);
imwrite(f.cdata, filename, 'Resolution', 300);

filename = sprintf('%s_index_check30.tiff', label );     
figureHandle = figure;
imagesc( ind_val( :,  :,   30 )<1.4);
f=getframe(figureHandle);
imwrite(f.cdata, filename, 'Resolution', 300);

filename = sprintf('%s_shell_check60.tiff', label );     
figureHandle = figure;
imagesc( dist3d(:,:,60)<0e-9&dist3d(:,:,60)>-2.1e-9 )
f=getframe(figureHandle);
imwrite(f.cdata, filename, 'Resolution', 300);

filename = sprintf('%s_shell_check30.tiff', label );     
figureHandle = figure;
imagesc( dist3d(:,:,30)<0e-9&dist3d(:,:,30)>-2.1e-9 )
f=getframe(figureHandle);
imwrite(f.cdata, filename, 'Resolution', 300);

filename = sprintf('%s_shell_check90.tiff', label );     
figureHandle = figure;
imagesc( dist3d(:,:,90)<0e-9&dist3d(:,:,90)>-2.1e-9 )
f=getframe(figureHandle);
imwrite(f.cdata, filename, 'Resolution', 300);


filename = sprintf('%s_dcheck60.tiff', label );     
figureHandle = figure;
imagesc( dist3d(:,:,60));
f=getframe(figureHandle);
imwrite(f.cdata, filename, 'Resolution', 300);
