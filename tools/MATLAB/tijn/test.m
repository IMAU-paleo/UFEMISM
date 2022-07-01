clc
clear all
close all

mesh = read_mesh_from_file('../../results_20220314_001/debug_ANT.nc');

d1 = ncread('../../../IMAU-ICE/results_20220314_001/debug_ANT.nc','dp_2D_01'  );
d2 = ncread('../../results_20220314_001/debug_ANT.nc'            ,'dp_2D_a_01');
H = plot_grid_and_mesh_data( d1,mesh,d2);
title(H.Ax1,'Hi');

d1 = ncread('../../../IMAU-ICE/results_20220314_001/debug_ANT.nc','dp_2D_02'  );
d2 = ncread('../../results_20220314_001/debug_ANT.nc'            ,'dp_2D_a_02');
H = plot_grid_and_mesh_data( d1,mesh,d2);
title(H.Ax1,'taudx');

d1 = ncread('../../../IMAU-ICE/results_20220314_001/debug_ANT.nc','dp_2D_03'  );
d2 = ncread('../../results_20220314_001/debug_ANT.nc'            ,'dp_2D_a_03');
H = plot_grid_and_mesh_data( d1,mesh,d2);
title(H.Ax1,'taudy');

d1 = ncread('../../../IMAU-ICE/results_20220314_001/debug_ANT.nc','dp_2D_04'  );
d2 = ncread('../../results_20220314_001/debug_ANT.nc'            ,'dp_2D_a_04');
H = plot_grid_and_mesh_data( d1,mesh,d2);
title(H.Ax1,'A');

d1 = ncread('../../../IMAU-ICE/results_20220314_001/debug_ANT.nc','dp_2D_05'  );
d2 = ncread('../../results_20220314_001/debug_ANT.nc'            ,'dp_2D_a_05');
H = plot_grid_and_mesh_data( d1,mesh,d2);
title(H.Ax1,'N');

d1 = ncread('../../../IMAU-ICE/results_20220314_001/debug_ANT.nc','dp_2D_07'  );
d2 = ncread('../../results_20220314_001/debug_ANT.nc'            ,'dp_2D_a_07');
d1 = min(d1,3e5);
d2 = min(d2,3e5);
H = plot_grid_and_mesh_data( d1,mesh,d2);
title(H.Ax1,'beta');

d1 = ncread('../../../IMAU-ICE/results_20220314_001/debug_ANT.nc','dp_2D_08'  );
d2 = ncread('../../results_20220314_001/debug_ANT.nc'            ,'dp_2D_a_08');
H = plot_grid_and_mesh_data( d1,mesh,d2);
title(H.Ax1,'F2');

d1 = ncread('../../../IMAU-ICE/results_20220314_001/debug_ANT.nc','dp_2D_09'  );
d2 = ncread('../../results_20220314_001/debug_ANT.nc'            ,'dp_2D_a_09');
H = plot_grid_and_mesh_data( d1,mesh,d2);
title(H.Ax1,'beta eff');

d1 = ncread('../../../IMAU-ICE/results_20220314_001/debug_ANT.nc','dp_2D_10'  );
d2 = ncread('../../results_20220314_001/debug_ANT.nc'            ,'dp_2D_a_10');
H = plot_grid_and_mesh_data( d1,mesh,d2);
title(H.Ax1,'u');

function H = plot_grid_and_mesh_data( d1,mesh,d2)

wa = 400;
ha = 400;

wf = 25 + wa + 25 + wa + 125;
hf = 35 + ha + 50;

dmin = min( min(d1(:)),min(d2(:)));
dmax = max( max(d1(:)),max(d2(:)));
if abs(dmin-1e-16)<1e-18 && abs(dmax-1e-16)<1e-18; dmin = 0.9e-16; dmax = 1.1e-16; end
if dmin == dmax; dmin = dmax-1; dmax = dmax+1; end
clim = [dmin,dmax];

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax1 = axes('parent',H.Fig,'units','pixels','position',[25      ,35,wa,ha],'xlim',[0,1],'ylim',[0,1],'fontsize',24,'clim',clim,'xtick',[],'ytick',[]);
H.Ax2 = axes('parent',H.Fig,'units','pixels','position',[25+wa+25,35,wa,ha],'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax],'clim',clim,...
  'fontsize',24,'clim',clim,'xtick',[],'ytick',[]);

H.Im1 = image('parent',H.Ax1,'xdata',[0,1],'ydata',[0,1],'cdata',d1','cdatamapping','scaled');
H.Im2 = patch('parent',H.Ax2,'vertices',mesh.V,'faces',mesh.Tri,'facevertexcdata',d2,'facecolor','interp','edgecolor','none');

pos = get(H.Ax2,'position');
H.Cbar = colorbar(H.Ax2,'location','eastoutside');
set(H.Ax2,'position',pos);

end