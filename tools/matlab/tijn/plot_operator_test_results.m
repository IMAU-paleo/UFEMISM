clc
clear all
close all

filename  = '../results_20211215_001/debug_ANT.nc';

mesh = ReadMeshFromFile( filename);

clim_map = [-1,1] * 1e-1;
clim_ddx = [-1,1] * 1e-6;

%% Read data

% a-grid (vertex)
e_aa   = ncread( filename, 'dp_2D_a_01');
e_ba   = ncread( filename, 'dp_2D_a_02');
e_ca   = ncread( filename, 'dp_2D_a_03');
edx_aa = ncread( filename, 'dp_2D_a_04');
edx_ba = ncread( filename, 'dp_2D_a_05');
edx_ca = ncread( filename, 'dp_2D_a_06');
edy_aa = ncread( filename, 'dp_2D_a_07');
edy_ba = ncread( filename, 'dp_2D_a_08');
edy_ca = ncread( filename, 'dp_2D_a_09');

% b-grid (triangle)
e_ab   = ncread( filename, 'dp_2D_b_01');
e_bb   = ncread( filename, 'dp_2D_b_02');
e_cb   = ncread( filename, 'dp_2D_b_03');
edx_ab = ncread( filename, 'dp_2D_b_04');
edx_bb = ncread( filename, 'dp_2D_b_05');
edx_cb = ncread( filename, 'dp_2D_b_06');
edy_ab = ncread( filename, 'dp_2D_b_07');
edy_bb = ncread( filename, 'dp_2D_b_08');
edy_cb = ncread( filename, 'dp_2D_b_09');

% c-grid (edge)
e_ac   = ncread( filename, 'dp_2D_c_01');
e_bc   = ncread( filename, 'dp_2D_c_02');
e_cc   = ncread( filename, 'dp_2D_c_03');
edx_ac = ncread( filename, 'dp_2D_c_04');
edx_bc = ncread( filename, 'dp_2D_c_05');
edx_cc = ncread( filename, 'dp_2D_c_06');
edy_ac = ncread( filename, 'dp_2D_c_07');
edy_bc = ncread( filename, 'dp_2D_c_08');
edy_cc = ncread( filename, 'dp_2D_c_09');

%% Mapping

% Plot
wa = 200;
ha = wa;

wf = 25 + 25 + wa + 25 + wa + 25 + wa + 25;
hf = wf + 25;

H1.Fig = figure('position',[200,200,wf,hf],'color','w','name','Mapping operators');

for i = 1:3
  for j = 1:3
    x = 25 + (j  )*25 + (j-1)*wa;
    y =      (4-i)*25 + (3-i)*ha;
    H1.Ax( i,j) = axes('parent',H1.Fig,'units','pixels','position',[x,y,wa,ha],...
      'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax],'xtick',[],'ytick',[],...
      'xaxislocation','top','fontsize',24,'clim',clim_map);
    H1.Axa( i,j) = axes('parent',H1.Fig,'units','pixels','position',[x,y,wa,ha],...
      'xtick',[],'ytick',[],'color','none','box','on');
  end
end

ylabel( H1.Ax(1,1),'from A');
ylabel( H1.Ax(2,1),'from B');
ylabel( H1.Ax(3,1),'from C');
xlabel( H1.Ax(1,1),'to A');
xlabel( H1.Ax(1,2),'to B');
xlabel( H1.Ax(1,3),'to C');

edgecolor = 'none';
ncols = 256;
cmap = parula(ncols);
markersize = 8;

patch('parent',H1.Ax( 1,1),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',e_aa,'facecolor','interp');
patch('parent',H1.Ax( 2,1),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',e_ba,'facecolor','interp');
patch('parent',H1.Ax( 3,1),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',e_ca,'facecolor','interp');

patch('parent',H1.Ax( 1,2),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',e_ab,'facecolor','flat');
patch('parent',H1.Ax( 2,2),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',e_bb,'facecolor','flat');
patch('parent',H1.Ax( 3,2),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',e_cb,'facecolor','flat');

d_c = e_ac;
ax  = H1.Ax( 1,3);
for ci = 1: ncols
  linedata(ci).x = [];
  linedata(ci).y = [];
end
for aci = 1: mesh.nAc
  clim = get(ax,'clim');
  ci = (d_c( aci) - clim(1)) / (clim(2) - clim(1));
  ci = max(1,min(ncols,1 + round(ci*(ncols-1))));
  linedata(ci).x( end+1) = mesh.VAc( aci,1);
  linedata(ci).y( end+1) = mesh.VAc( aci,2);
end
for ci = 1: ncols
  line('parent',ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
    'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',markersize);
end

d_c = e_bc;
ax  = H1.Ax( 2,3);
for ci = 1: ncols
  linedata(ci).x = [];
  linedata(ci).y = [];
end
for aci = 1: mesh.nAc
  clim = get(ax,'clim');
  ci = (d_c( aci) - clim(1)) / (clim(2) - clim(1));
  ci = max(1,min(ncols,1 + round(ci*(ncols-1))));
  linedata(ci).x( end+1) = mesh.VAc( aci,1);
  linedata(ci).y( end+1) = mesh.VAc( aci,2);
end
for ci = 1: ncols
  line('parent',ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
    'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',markersize);
end

d_c = e_cc;
ax  = H1.Ax( 3,3);
for ci = 1: ncols
  linedata(ci).x = [];
  linedata(ci).y = [];
end
for aci = 1: mesh.nAc
  clim = get(ax,'clim');
  ci = (d_c( aci) - clim(1)) / (clim(2) - clim(1));
  ci = max(1,min(ncols,1 + round(ci*(ncols-1))));
  linedata(ci).x( end+1) = mesh.VAc( aci,1);
  linedata(ci).y( end+1) = mesh.VAc( aci,2);
end
for ci = 1: ncols
  line('parent',ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
    'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',markersize);
end

%% d/dx

% Plot
wa = 200;
ha = wa;

wf = 25 + 25 + wa + 25 + wa + 25 + wa + 25;
hf = wf + 25;

H2.Fig = figure('position',[400,200,wf,hf],'color','w','name','d/dx operators');

for i = 1:3
  for j = 1:3
    x = 25 + (j  )*25 + (j-1)*wa;
    y =      (4-i)*25 + (3-i)*ha;
    H2.Ax( i,j) = axes('parent',H2.Fig,'units','pixels','position',[x,y,wa,ha],...
      'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax],'xtick',[],'ytick',[],...
      'xaxislocation','top','fontsize',24,'clim',clim_ddx);
    H2.Axa( i,j) = axes('parent',H2.Fig,'units','pixels','position',[x,y,wa,ha],...
      'xtick',[],'ytick',[],'color','none','box','on');
  end
end

ylabel( H2.Ax(1,1),'from A');
ylabel( H2.Ax(2,1),'from B');
ylabel( H2.Ax(3,1),'from C');
xlabel( H2.Ax(1,1),'to A');
xlabel( H2.Ax(1,2),'to B');
xlabel( H2.Ax(1,3),'to C');

edgecolor = 'none';
ncols = 256;
cmap = parula(ncols);
markersize = 8;

patch('parent',H2.Ax( 1,1),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',edx_aa,'facecolor','interp');
patch('parent',H2.Ax( 2,1),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',edx_ba,'facecolor','interp');
patch('parent',H2.Ax( 3,1),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',edx_ca,'facecolor','interp');

patch('parent',H2.Ax( 1,2),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',edx_ab,'facecolor','flat');
patch('parent',H2.Ax( 2,2),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',edx_bb,'facecolor','flat');
patch('parent',H2.Ax( 3,2),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',edx_cb,'facecolor','flat');

d_c = edx_ac;
ax  = H2.Ax( 1,3);
for ci = 1: ncols
  linedata(ci).x = [];
  linedata(ci).y = [];
end
for aci = 1: mesh.nAc
  clim = get(ax,'clim');
  ci = (d_c( aci) - clim(1)) / (clim(2) - clim(1));
  ci = max(1,min(ncols,1 + round(ci*(ncols-1))));
  linedata(ci).x( end+1) = mesh.VAc( aci,1);
  linedata(ci).y( end+1) = mesh.VAc( aci,2);
end
for ci = 1: ncols
  line('parent',ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
    'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',markersize);
end

d_c = edx_bc;
ax  = H2.Ax( 2,3);
for ci = 1: ncols
  linedata(ci).x = [];
  linedata(ci).y = [];
end
for aci = 1: mesh.nAc
  clim = get(ax,'clim');
  ci = (d_c( aci) - clim(1)) / (clim(2) - clim(1));
  ci = max(1,min(ncols,1 + round(ci*(ncols-1))));
  linedata(ci).x( end+1) = mesh.VAc( aci,1);
  linedata(ci).y( end+1) = mesh.VAc( aci,2);
end
for ci = 1: ncols
  line('parent',ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
    'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',markersize);
end

d_c = edx_cc;
ax  = H2.Ax( 3,3);
for ci = 1: ncols
  linedata(ci).x = [];
  linedata(ci).y = [];
end
for aci = 1: mesh.nAc
  clim = get(ax,'clim');
  ci = (d_c( aci) - clim(1)) / (clim(2) - clim(1));
  ci = max(1,min(ncols,1 + round(ci*(ncols-1))));
  linedata(ci).x( end+1) = mesh.VAc( aci,1);
  linedata(ci).y( end+1) = mesh.VAc( aci,2);
end
for ci = 1: ncols
  line('parent',ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
    'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',markersize);
end

%% d/dy

% Plot
wa = 200;
ha = wa;

wf = 25 + 25 + wa + 25 + wa + 25 + wa + 25;
hf = wf + 25;

H3.Fig = figure('position',[600,200,wf,hf],'color','w','name','d/dy operators');

for i = 1:3
  for j = 1:3
    x = 25 + (j  )*25 + (j-1)*wa;
    y =      (4-i)*25 + (3-i)*ha;
    H3.Ax( i,j) = axes('parent',H3.Fig,'units','pixels','position',[x,y,wa,ha],...
      'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax],'xtick',[],'ytick',[],...
      'xaxislocation','top','fontsize',24,'clim',clim_ddx);
    H3.Axa( i,j) = axes('parent',H3.Fig,'units','pixels','position',[x,y,wa,ha],...
      'xtick',[],'ytick',[],'color','none','box','on');
  end
end

ylabel( H3.Ax(1,1),'from A');
ylabel( H3.Ax(2,1),'from B');
ylabel( H3.Ax(3,1),'from C');
xlabel( H3.Ax(1,1),'to A');
xlabel( H3.Ax(1,2),'to B');
xlabel( H3.Ax(1,3),'to C');

edgecolor = 'none';
ncols = 256;
cmap = parula(ncols);
markersize = 8;

patch('parent',H3.Ax( 1,1),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',edy_aa,'facecolor','interp');
patch('parent',H3.Ax( 2,1),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',edy_ba,'facecolor','interp');
patch('parent',H3.Ax( 3,1),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',edy_ca,'facecolor','interp');

patch('parent',H3.Ax( 1,2),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',edy_ab,'facecolor','flat');
patch('parent',H3.Ax( 2,2),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',edy_bb,'facecolor','flat');
patch('parent',H3.Ax( 3,2),'vertices',mesh.V,'faces',mesh.Tri,'edgecolor',edgecolor,...
  'facevertexcdata',edy_cb,'facecolor','flat');

d_c = edy_ac;
ax  = H3.Ax( 1,3);
for ci = 1: ncols
  linedata(ci).x = [];
  linedata(ci).y = [];
end
for aci = 1: mesh.nAc
  clim = get(ax,'clim');
  ci = (d_c( aci) - clim(1)) / (clim(2) - clim(1));
  ci = max(1,min(ncols,1 + round(ci*(ncols-1))));
  linedata(ci).x( end+1) = mesh.VAc( aci,1);
  linedata(ci).y( end+1) = mesh.VAc( aci,2);
end
for ci = 1: ncols
  line('parent',ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
    'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',markersize);
end

d_c = edy_bc;
ax  = H3.Ax( 2,3);
for ci = 1: ncols
  linedata(ci).x = [];
  linedata(ci).y = [];
end
for aci = 1: mesh.nAc
  clim = get(ax,'clim');
  ci = (d_c( aci) - clim(1)) / (clim(2) - clim(1));
  ci = max(1,min(ncols,1 + round(ci*(ncols-1))));
  linedata(ci).x( end+1) = mesh.VAc( aci,1);
  linedata(ci).y( end+1) = mesh.VAc( aci,2);
end
for ci = 1: ncols
  line('parent',ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
    'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',markersize);
end

d_c = edy_cc;
ax  = H3.Ax( 3,3);
for ci = 1: ncols
  linedata(ci).x = [];
  linedata(ci).y = [];
end
for aci = 1: mesh.nAc
  clim = get(ax,'clim');
  ci = (d_c( aci) - clim(1)) / (clim(2) - clim(1));
  ci = max(1,min(ncols,1 + round(ci*(ncols-1))));
  linedata(ci).x( end+1) = mesh.VAc( aci,1);
  linedata(ci).y( end+1) = mesh.VAc( aci,2);
end
for ci = 1: ncols
  line('parent',ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
    'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',markersize);
end