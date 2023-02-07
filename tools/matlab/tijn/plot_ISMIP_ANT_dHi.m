clc
clear all
close all

foldername = '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_1950-1995_oceanccsm4_rcp8.5';
% foldername = '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_1995-2100_oceanccsm4_rcp8.5';

filename = [foldername '/restart_ANT_00001.nc'];

mesh = read_mesh_from_file( filename);

clim = [-400,400];

%% Set up GUI

wa = 800;
ha = 800;

margin_left   = 25;
margin_right  = 100;
margin_bottom = 25;
margin_top    = 50;

wf = margin_left   + wa + margin_right;
hf = margin_bottom + ha + margin_top;

H.Fig = figure('position',[200,100,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margin_left,margin_bottom,wa,ha],...
  'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax],'fontsize',24,...
  'xtick',[],'ytick',[],'box','on','clim',clim);
H.dHi = patch('parent',H.Ax,'vertices',mesh.V,'faces',mesh.Tri,'facevertexcdata',zeros(mesh.nV,1),...
  'facecolor','interp','edgecolor','none');

cmap = bluewhiteredmap(255);
colormap(H.Ax,cmap);
H.Cbar = colorbar( H.Ax,'location','eastoutside');
ylabel(H.Cbar,'\Delta H (m)')

%% Plot PD grounding line and ice front
load('ANT_GL_hires.mat')
load('ANT_CF_hires.mat')
line('parent',H.Ax,'xdata',CF.x,'ydata',CF.y);
line('parent',H.Ax,'xdata',GL.x,'ydata',GL.y);

%% Make animation
Hi_init = ncread( filename,'Hi',[1,1],[Inf,1]);
time = ncread( filename,'time');
for ti = 1: length( time)
  title(H.Ax,['Year: ' num2str( time(ti))])
  Hi = ncread( filename,'Hi',[1,ti],[Inf,1]);
  dHi = Hi - Hi_init;
  set(H.dHi,'facevertexcdata',dHi);
  drawnow('update');
  pause(0.1);
end