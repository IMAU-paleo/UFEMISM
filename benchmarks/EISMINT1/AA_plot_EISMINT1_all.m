clc
clear all
close all

foldernames = {
  'EISMINT1_A_SIASSA';
  'EISMINT1_B_SIASSA';
  'EISMINT1_C_SIASSA';
  'EISMINT1_D_SIASSA';
  'EISMINT1_E_SIASSA';
  'EISMINT1_F_SIASSA';
  };

% Plot parameters
wa = 800;
ha = 500;

% Read UFEMISM output
results = read_data( foldernames);
% save('tempdata.mat','results')
% load('tempdata.mat')

% Plot all figures
plot_Fig_2a( 'Huybrechts1996_figures_raw/Fig_2a.png',results,wa,ha)
plot_Fig_2c( 'Huybrechts1996_figures_raw/Fig_2c.png',results,wa,ha)
plot_Fig_2e( 'Huybrechts1996_figures_raw/Fig_2e.png',results,wa,ha)
plot_Fig_2f( 'Huybrechts1996_figures_raw/Fig_2f.png',results,wa,ha)
plot_Fig_3a( 'Huybrechts1996_figures_raw/Fig_3a.png',results,wa,ha)
plot_Fig_3b( 'Huybrechts1996_figures_raw/Fig_3b.png',results,wa,ha)
plot_Fig_3c( 'Huybrechts1996_figures_raw/Fig_3c.png',results,wa,ha)
plot_Fig_3d( 'Huybrechts1996_figures_raw/Fig_3d.png',results,wa,ha)
plot_Fig_4a( 'Huybrechts1996_figures_raw/Fig_4a.png',results,wa,ha)
plot_Fig_4c( 'Huybrechts1996_figures_raw/Fig_4c.png',results,wa,ha)
plot_Fig_4d( 'Huybrechts1996_figures_raw/Fig_4d.png',results,wa,ha)

plot_Fig_5a( 'Huybrechts1996_figures_raw/Fig_5a.png',results,wa,ha)
plot_Fig_5b( 'Huybrechts1996_figures_raw/Fig_5b.png',results,wa,ha)
plot_Fig_5c( 'Huybrechts1996_figures_raw/Fig_5c.png',results,wa,ha)
plot_Fig_5d( 'Huybrechts1996_figures_raw/Fig_5d.png',results,wa,ha)
plot_Fig_6a( 'Huybrechts1996_figures_raw/Fig_6a.png',results,wa,ha)
plot_Fig_6b( 'Huybrechts1996_figures_raw/Fig_6b.png',results,wa,ha)
plot_Fig_6c( 'Huybrechts1996_figures_raw/Fig_6c.png',results,wa,ha)
plot_Fig_6d( 'Huybrechts1996_figures_raw/Fig_6d.png',results,wa,ha)
plot_Fig_7a( 'Huybrechts1996_figures_raw/Fig_7a.png',results,wa,ha)
plot_Fig_7b( 'Huybrechts1996_figures_raw/Fig_7b.png',results,wa,ha)
plot_Fig_7c( 'Huybrechts1996_figures_raw/Fig_7c.png',results,wa,ha)
plot_Fig_7d( 'Huybrechts1996_figures_raw/Fig_7d.png',results,wa,ha)

function results = read_data( foldernames)

for fi = 1: length( foldernames)
  
  filename_restart     = [foldernames{ fi} '/restart_ANT_00001.nc'    ];
  filename_help_fields = [foldernames{ fi} '/help_fields_ANT_00001.nc'];
  
  results( fi).name = foldernames{ fi}(1:10);
  
  results( fi).time = ncread( filename_restart,'time');
  
  zeta = ncread( filename_restart,'zeta');
  nz = length( zeta);
  
  results( fi).zeta = zeta;
  results( fi).nz   = nz;
  
  % Read mesh
  mesh = read_mesh_from_file( filename_restart);
  
  %% Transect
  
  % Calculate transect matrix
  trans.n = 100;
  trans.x = linspace( 0, mesh.xmax, trans.n)';
  trans.y = zeros( size( trans.x));
  trans.A = calc_transect_matrix_a( mesh, trans.x, trans.y);
  
  % Read model output
  time = ncread( filename_restart,'time');
  ti = length( time);
  
  Hi        = ncread( filename_restart    , 'Hi'    ,[1  ,ti],[Inf    ,1]);
  u_3D      = ncread( filename_help_fields, 'u_3D'  ,[1,1,ti],[Inf,Inf,1]);
  Ti_3D     = ncread( filename_restart    , 'Ti'    ,[1,1,ti],[Inf,Inf,1]);
  Ti_pmp_3D = ncread( filename_help_fields, 'Ti_pmp',[1,1,ti],[Inf,Inf,1]);
  Ti_hom_3D = Ti_3D - Ti_pmp_3D;
  for ti = 1: size( Ti_hom_3D,3)
    m = Hi( :,ti) == 0;
    for k = 1: size( Ti_hom_3D,2)
      henk = Ti_hom_3D( :,k,ti);
      henk( m) = 0;
      Ti_hom_3D( :,k,ti) = henk;
    end
  end
  
  % Calculate transects of ice thickness, velocity, and temperature
  trans.Hi = trans.A * Hi;
  
  trans.u_3D      = zeros( trans.n,nz);
  trans.Ti_hom_3D = zeros( trans.n,nz);
  for k = 1: nz
    trans.u_3D(      :,k) = trans.A * u_3D(      :,k);
    trans.Ti_hom_3D( :,k) = trans.A * Ti_hom_3D( :,k);
  end
  
  % Vertically averaged velocity and ice flux
  trans.u_vav = zeros( trans.n,1);
  trans.Q_ice = zeros( trans.n,1);
  for i = 1: trans.n
    trans.u_vav( i) = vertical_average( zeta, trans.u_3D( i,:));
    trans.Q_ice( i) = trans.u_vav( i) * trans.Hi( i);
  end
  
  results( fi).trans = trans;
  
  %% Margin position
  results( fi).R_margin = zeros( size( time));
  for ti = 1: length( time)
    Hi  = ncread( filename_restart, 'Hi',[1,ti],[Inf,1]);
    n = 0;
    S = 0;
    for vi = 1: mesh.nV
      for ci = 1: mesh.nC( vi)
        vj = mesh.C( vi,ci);
        if ((Hi( vi) == 0 && Hi( vj) > 0) || (Hi( vi) > 0 && Hi( vj) == 0))
          r = norm( (mesh.V( vi,:) + mesh.V( vj,:))/2);
          n = n + 1;
          S = S + r;
        end
      end
    end
    results( fi).R_margin( ti) = S / n;
  end
 
  %% Summit
  
  results( fi).summit.time = time;
  
  % Find triangle containing the midpoint
  p = [0,0];
  mesh.TriMap    = zeros( mesh.nTri,1);
  mesh.TriStack1 = zeros( mesh.nTri,1);
  mesh.TriStack2 = zeros( mesh.nTri,1);
  ti = find_containing_triangle( mesh, p, 1);
  
  % The three vertices spanning the triangle
  via = mesh.Tri( ti,1);
  vib = mesh.Tri( ti,2);
  vic = mesh.Tri( ti,3);
  pa  = mesh.V( via,:);
  pb  = mesh.V( vib,:);
  pc  = mesh.V( vic,:);
  
  % The three interpolation weights
  Atri_abp = calc_triangle_area( pa, pb, p);
  Atri_bcp = calc_triangle_area( pb, pc, p);
  Atri_cap = calc_triangle_area( pc, pa, p);
  Atri_abc = Atri_abp + Atri_bcp + Atri_cap;

  wa = Atri_bcp / Atri_abc;
  wb = Atri_cap / Atri_abc;
  wc = Atri_abp / Atri_abc;
  
  % Read data
  Hi        = ncread( filename_restart    , 'Hi'    );
  w_3D      = ncread( filename_help_fields, 'w_3D'  );
  Ti_3D     = ncread( filename_restart    , 'Ti'    );
  Ti_pmp_3D = ncread( filename_help_fields, 'Ti_pmp');
  Ti_hom_3D = Ti_3D - Ti_pmp_3D;
  
  results( fi).summit.Hi        =         (wa * Hi(        via,:  ) + wb * Hi(        vib,:  ) + wc * Hi(        vic,:  ))';
  results( fi).summit.w_3D      = permute( wa * w_3D(      via,:,:) + wb * w_3D(      vib,:,:) + wc * w_3D(      vic,:,:), [2,3,1]);
  results( fi).summit.Ti_3D     = permute( wa * Ti_3D(     via,:,:) + wb * Ti_3D(     vib,:,:) + wc * Ti_3D(     vic,:,:), [2,3,1]);
  results( fi).summit.Ti_hom_3D = permute( wa * Ti_hom_3D( via,:,:) + wb * Ti_hom_3D( vib,:,:) + wc * Ti_hom_3D( vic,:,:), [2,3,1]);
 
  %% Midpoint
  
  results( fi).midpoint.time = time;
  
  % Find triangle containing the midpoint
  p = [mesh.xmax / 2 + 25e3,0];
  mesh.TriMap    = zeros( mesh.nTri,1);
  mesh.TriStack1 = zeros( mesh.nTri,1);
  mesh.TriStack2 = zeros( mesh.nTri,1);
  ti = find_containing_triangle( mesh, p, 1);
  
  % The three vertices spanning the triangle
  via = mesh.Tri( ti,1);
  vib = mesh.Tri( ti,2);
  vic = mesh.Tri( ti,3);
  pa  = mesh.V( via,:);
  pb  = mesh.V( vib,:);
  pc  = mesh.V( vic,:);
  
  % The three interpolation weights
  Atri_abp = calc_triangle_area( pa, pb, p);
  Atri_bcp = calc_triangle_area( pb, pc, p);
  Atri_cap = calc_triangle_area( pc, pa, p);
  Atri_abc = Atri_abp + Atri_bcp + Atri_cap;

  wa = Atri_bcp / Atri_abc;
  wb = Atri_cap / Atri_abc;
  wc = Atri_abp / Atri_abc;
  
  % Read data
  Hi        = ncread( filename_restart    , 'Hi'    );
  u_3D      = ncread( filename_help_fields, 'u_3D'  );
  w_3D      = ncread( filename_help_fields, 'w_3D'  );
  Ti_3D     = ncread( filename_restart    , 'Ti'    );
  Ti_pmp_3D = ncread( filename_help_fields, 'Ti_pmp');
  Ti_hom_3D = Ti_3D - Ti_pmp_3D;
  
  results( fi).midpoint.Hi        =         (wa * Hi(        via,:  ) + wb * Hi(        vib,:  ) + wc * Hi(        vic,:  ))';
  results( fi).midpoint.u_3D      = permute( wa * u_3D(      via,:,:) + wb * u_3D(      vib,:,:) + wc * u_3D(      vic,:,:), [2,3,1]);
  results( fi).midpoint.w_3D      = permute( wa * w_3D(      via,:,:) + wb * w_3D(      vib,:,:) + wc * w_3D(      vic,:,:), [2,3,1]);
  results( fi).midpoint.Ti_hom_3D = permute( wa * Ti_hom_3D( via,:,:) + wb * Ti_hom_3D( vib,:,:) + wc * Ti_hom_3D( vic,:,:), [2,3,1]);
  
  results( fi).midpoint.u_vav  = zeros( size( results( fi).midpoint.Hi));
  for ti = 1: length( results( fi).midpoint.time)
    results( fi).midpoint.u_vav( ti) = vertical_average( zeta, results( fi).midpoint.u_3D( :,ti));
  end
  results( fi).midpoint.Q_ice = results( fi).midpoint.Hi .* results( fi).midpoint.u_vav;
  
end

end

function plot_Fig_2a( figname, results, wa, ha)
% Plot ice thickness transect for the fixed-margin steady-state experiment EISMINT1_D

%% Set up GUI

xlim = [-120,765];
ylim = [-630,4200];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 2a');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:100:700
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:1000:4000
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_D')
    R = results( fi);
    line('xdata',R.trans.x / 1e3,'ydata',R.trans.Hi,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_2c( figname, results, wa, ha)
% Plot ice flux transect for the fixed-margin steady-state experiment EISMINT1_D

%% Set up GUI

xlim = [-170,755];
ylim = [-60,370] * 1e3;

figure('position',[200,200,wa,ha],'color','w','name','Fig. 2c');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:100:700
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:50:350
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_D')
    R = results( fi);
    line('xdata',R.trans.x / 1e3,'ydata',R.trans.Q_ice,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_2e( figname, results, wa, ha)
% Plot vertically averaged velocity transect for the fixed-margin steady-state experiment EISMINT1_D

%% Set up GUI

xlim = [-105,755];
ylim = [-25,160];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 2e');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:100:700
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:20:150
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_D')
    R = results( fi);
    line('xdata',R.trans.x / 1e3,'ydata',R.trans.u_vav,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_2f( figname, results, wa, ha)
% Plot vertical profile of horizontal velocity at the midpoint
% for the fixed-margin steady-state experiment EISMINT1_D

%% Set up GUI

xlim = [-5.2,41];
ylim = [-0.17,1.06];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 2f');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:10:40
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:0.2:1
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% % Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_D')
    R = results( fi);
    line('xdata',R.midpoint.u_3D(:,end),'ydata',1-R.zeta,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_3a( figname, results, wa, ha)
% Plot homologous basal temperature transect for the fixed-margin steady-state experiment EISMINT1_D

%% Set up GUI

xlim = [-120,775];
ylim = [-13.2,2.5];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 3a');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:100:700
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = -10:2:2
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_D')
    R = results( fi);
    line('xdata',R.trans.x / 1e3,'ydata',R.trans.Ti_hom_3D( :,end),'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_3b( figname, results, wa, ha)
% Plot vertical profile of vertical velocity at the ice divide
% for the fixed-margin steady-state experiment EISMINT1_D

%% Set up GUI

xlim = [-0.46,0.11];
ylim = [-0.17,1.06];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 3b');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = -0.4:0.1:0.1
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:0.2:1
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_D')
    R = results( fi);
    line('xdata',R.summit.w_3D(:,end),'ydata',1-R.zeta,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_3c( figname, results, wa, ha)
% Plot vertical profile of homologous temperature at the ice divide
% for the fixed-margin steady-state experiment EISMINT1_D

%% Set up GUI

xlim = [-39.8,0.6];
ylim = [-0.17,1.06];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 3c');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = -35:5:0
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:0.2:1
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_D')
    R = results( fi);
    line('xdata',R.summit.Ti_hom_3D(:,end),'ydata',1-R.zeta,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_3d( figname, results, wa, ha)
% Plot vertical profile of homologous temperature at the midpoint
% for the fixed-margin steady-state experiment EISMINT1_D

%% Set up GUI

xlim = [-39.8,1];
ylim = [-0.17,1.04];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 3d');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = -35:5:0
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:0.2:1
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_D')
    R = results( fi);
    line('xdata',R.midpoint.Ti_hom_3D(:,end),'ydata',1-R.zeta,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_4a( figname, results, wa, ha)
% Plot time series of ice thickness change at the divide
% for the fixed-margin glacial cycle experiments (EISMINT1_E and EISMINT1_F)

%% Set up GUI

xlim = [-35,212] * 1e3;
ylim = [-740,345];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 4a');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:50*1e3:200*1e3
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = -600:200:300
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_E')
    R = results( fi);
    line('xdata',R.time,'ydata',R.summit.Hi - R.summit.Hi(1),'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
  if strcmpi( results( fi).name,'EISMINT1_F')
    R = results( fi);
    line('xdata',R.time,'ydata',R.summit.Hi - R.summit.Hi(1),'linestyle','none','color','b','marker','o','markerfacecolor','b');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_4c( figname, results, wa, ha)
% Plot time series of ice flux at the midpoint
% for the fixed-margin glacial cycle experiments (EISMINT1_E and EISMINT1_F)

%% Set up GUI

xlim = [-53,213] * 1e3;
ylim = [-158,235] * 1e3;

figure('position',[200,200,wa,ha],'color','w','name','Fig. 4c');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:50*1e3:200*1e3
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = -100*1e3:50*1e3:200*1e3
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_E')
    R = results( fi);
    line('xdata',R.time,'ydata',R.midpoint.Q_ice,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
  if strcmpi( results( fi).name,'EISMINT1_F')
    R = results( fi);
    line('xdata',R.time,'ydata',R.midpoint.Q_ice,'linestyle','none','color','b','marker','o','markerfacecolor','b');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_4d( figname, results, wa, ha)
% Plot time series of basal temperature change at the ice divide
% for the fixed-margin glacial cycle experiments (EISMINT1_E and EISMINT1_F)

%% Set up GUI

xlim = [-28,212] * 1e3;
ylim = [-2.6,5.35];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 4d');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:50*1e3:200*1e3
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = -1:5
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_E')
    R = results( fi);
    line('xdata',R.time,'ydata',R.summit.Ti_3D(end,:) - R.summit.Ti_3D(end,1),'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
  if strcmpi( results( fi).name,'EISMINT1_F')
    R = results( fi);
    line('xdata',R.time,'ydata',R.summit.Ti_3D(end,:) - R.summit.Ti_3D(end,1),'linestyle','none','color','b','marker','o','markerfacecolor','b');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end

function plot_Fig_5a( figname, results, wa, ha)
% Plot ice thickness transect for the moving-margin steady-state experiment EISMINT1_A

%% Set up GUI

xlim = [-125,765];
ylim = [-570,3250];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 5a');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:100:700
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:500:3000
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_A')
    R = results( fi);
    line('xdata',R.trans.x / 1e3,'ydata',R.trans.Hi,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_5b( figname, results, wa, ha)
% Plot ice flux transect for the moving-margin steady-state experiment EISMINT1_A

%% Set up GUI

xlim = [-177,770];
ylim = [-25,150] * 1e3;

figure('position',[200,200,wa,ha],'color','w','name','Fig. 5b');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:100:700
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:20*1e3:140*1e3
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_A')
    R = results( fi);
    line('xdata',R.trans.x / 1e3,'ydata',R.trans.Q_ice,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_5c( figname, results, wa, ha)
% Plot ice velocity transect for the moving-margin steady-state experiment EISMINT1_A

%% Set up GUI

xlim = [-90,760];
ylim = [-15,85];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 5c');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:100:700
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:20:80
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_A')
    R = results( fi);
    line('xdata',R.trans.x / 1e3,'ydata',R.trans.u_vav,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_5d( figname, results, wa, ha)
% Plot vertical profile of horizontal velocity at the midpoint
% for the moving-margin steady-state experiment EISMINT1_A

%% Set up GUI

xlim = [-9.5,71.5];
ylim = [-0.17,1.06];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 5d');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:10:70
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:0.2:1
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_A')
    R = results( fi);
    line('xdata',R.midpoint.u_3D(:,end),'ydata',1-R.zeta,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_6a( figname, results, wa, ha)
% Plot homologous basal temperature transect for the moving-margin steady-state experiment EISMINT1_A

%% Set up GUI

xlim = [-120,775];
ylim = [-18,2.7];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 6a');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:100:700
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = -14:2:2
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_A')
    R = results( fi);
    line('xdata',R.trans.x / 1e3,'ydata',R.trans.Ti_hom_3D( :,end),'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_6b( figname, results, wa, ha)
% Plot vertical profile of vertical velocity at the ice divide
% for the moving-margin steady-state experiment EISMINT1_A

%% Set up GUI

xlim = [-0.695,0.12];
ylim = [-0.17,1.06];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 6b');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = -0.6:0.1:0.1
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:0.2:1
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_A')
    R = results( fi);
    line('xdata',R.summit.w_3D(:,end),'ydata',1-R.zeta,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_6c( figname, results, wa, ha)
% Plot vertical profile of homologous temperature at the ice divide
% for the moving-margin steady-state experiment EISMINT1_A

%% Set up GUI

xlim = [-39.8,1];
ylim = [-0.17,1.06];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 6c');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = -35:5:0
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:0.2:1
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_A')
    R = results( fi);
    line('xdata',R.summit.Ti_hom_3D(:,end),'ydata',1-R.zeta,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_6d( figname, results, wa, ha)
% Plot vertical profile of homologous temperature at the midpoint
% for the moving-margin steady-state experiment EISMINT1_A

%% Set up GUI

xlim = [-34.5,6];
ylim = [-0.17,1.06];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 6d');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = -30:5:5
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 0:0.2:1
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_A')
    R = results( fi);
    line('xdata',R.midpoint.Ti_hom_3D(:,end),'ydata',1-R.zeta,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_7a( figname, results, wa, ha)
% Plot time series of ice thickness change at the divide
% for the moving-margin glacial cycle experiments (EISMINT1_B and EISMINT1_C)

%% Set up GUI

xlim = [-33,213] * 1e3;
ylim = [-740,395];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 7a');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:50*1e3:200*1e3
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = -600:200:200
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_B')
    R = results( fi);
    line('xdata',R.time,'ydata',R.summit.Hi - R.summit.Hi(1),'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
  if strcmpi( results( fi).name,'EISMINT1_C')
    R = results( fi);
    line('xdata',R.time,'ydata',R.summit.Hi - R.summit.Hi(1),'linestyle','none','color','b','marker','o','markerfacecolor','b');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_7b( figname, results, wa, ha)
% Plot time series of mass flux at midpoint
% for the moving-margin glacial cycle experiments (EISMINT1_B and EISMINT1_C)

%% Set up GUI

xlim = [-45,215] * 1e3;
ylim = [25,134] * 1e3;

figure('position',[200,200,wa,ha],'color','w','name','Fig. 7b');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:50*1e3:200*1e3
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 40*1e3:20*1e3:120*1e3
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_B')
    R = results( fi);
    line('xdata',R.time,'ydata',R.midpoint.Q_ice,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
  if strcmpi( results( fi).name,'EISMINT1_C')
    R = results( fi);
    line('xdata',R.time,'ydata',R.midpoint.Q_ice,'linestyle','none','color','b','marker','o','markerfacecolor','b');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_7c( figname, results, wa, ha)
% Plot time series of basal temperature change at the ice divide
% for the moving-margin glacial cycle experiments (EISMINT1_B and EISMINT1_C)

%% Set up GUI

xlim = [-26,213] * 1e3;
ylim = [-9,5.8];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 7c');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:50*1e3:200*1e3
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = -6:2:4
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_B')
    R = results( fi);
    line('xdata',R.time,'ydata',R.summit.Ti_3D(end,:) - R.summit.Ti_3D(end,1),'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
  if strcmpi( results( fi).name,'EISMINT1_C')
    R = results( fi);
    line('xdata',R.time,'ydata',R.summit.Ti_3D(end,:) - R.summit.Ti_3D(end,1),'linestyle','none','color','b','marker','o','markerfacecolor','b');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end
function plot_Fig_7d( figname, results, wa, ha)
% Plot time series of margin position
% for the moving-margin glacial cycle experiments (EISMINT1_B and EISMINT1_C)

%% Set up GUI

xlim = [-26,213] * 1e3;
ylim = [22.7,31.3];

figure('position',[200,200,wa,ha],'color','w','name','Fig. 7d');
axes('position',[0,0,1,1],'xtick',[],'ytick',[],'xlim',xlim,'ylim',ylim);
im = imread( figname);
image('cdata',im,'xdata',xlim,'ydata',fliplr(ylim));

% % Grid
% for x = 0:50*1e3:200*1e3
%   line('xdata',[0,0]+x,'ydata',ylim,'linestyle','-','linewidth',2,'color','k')
% end
% for y = 24:31
%   line('ydata',[0,0]+y,'xdata',xlim,'linestyle','-','linewidth',2,'color','k')
% end

% Plot UFEMISM data
for fi = 1: length( results)
  if strcmpi( results( fi).name,'EISMINT1_B')
    R = results( fi);
    R_margin = R.R_margin / 50000 + 16;
    line('xdata',R.time,'ydata',R_margin,'linestyle','none','color','r','marker','o','markerfacecolor','r');
  end
  if strcmpi( results( fi).name,'EISMINT1_C')
    R = results( fi);
    R_margin = R.R_margin / 50000 + 16;
    line('xdata',R.time,'ydata',R_margin,'linestyle','none','color','b','marker','o','markerfacecolor','b');
  end
end

% Save to file
i = strfind( figname,'/');
filename = ['AA_' figname(i+1:end)];
if exist(filename,'file')
  delete(filename)
end
saveas(gcf,filename);

end