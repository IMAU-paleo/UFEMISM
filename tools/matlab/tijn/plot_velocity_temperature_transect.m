clc
clear all
close all

foldername = '../../../benchmarks/EISMINT1/EISMINT1_A_SIASSA';

xlim     = [-700,700];
ylim     = [0,3500];
clim     = [220,275];
trans.n  = 60;
dt_pause = 0.1;

%% Read and process data

% Read mesh
filename_restart     = [foldername '/restart_ANT_00001.nc'];
filename_help_fields = [foldername '/help_fields_ANT_00001.nc'];
mesh = read_mesh_from_file( filename_help_fields);

% Calculate transect matrix
trans.x = linspace( mesh.xmin, mesh.xmax, trans.n)';
trans.y = zeros( size( trans.x));
trans.A = calc_transect_matrix_a( mesh, trans.x, trans.y);

zeta = ncread( filename_help_fields, 'zeta');
nz = length( zeta);

%% Set up GUI

wa = 800;
ha = 500;
margins_hor = [150,25];
margins_ver = [80 ,50];

wf = margins_hor( 1) + wa + margins_hor( 2);
hf = margins_ver( 1) + ha + margins_ver( 2);

H.Fig = figure( 'color','w','position',[200,200,wf,hf]);
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margins_hor(1),margins_ver(1),wa,ha],'fontsize',24,...
  'xlim',xlim,'ylim',ylim,'clim',clim);
xlabel(H.Ax,'x (km)');
ylabel(H.Ax,'z (m)');

H.Patch   = patch('parent',H.Ax,'vertices',[],'faces',[],'facevertexcdata',[],'facecolor','interp','edgecolor','none');
H.IceSurf = line( 'parent',H.Ax,'xdata',[],'ydata',[],'linewidth',3);

colormap( H.Ax,parula(16))

%% Loop over time
time = ncread( filename_help_fields, 'time');
Quiv = [];

for ti = length( time)
  
  title(H.Ax,['Time = ' num2str( time(ti)) ' yr']);

  %% Calculate H,Ti,u,v,w transects
  Hi   = ncread( filename_restart    ,'Hi'  ,[1,  ti],[Inf,    1]);
  Ti   = ncread( filename_restart    ,'Ti'  ,[1,1,ti],[Inf,Inf,1]);
  u_3D = ncread( filename_help_fields,'u_3D',[1,1,ti],[Inf,Inf,1]);
  v_3D = ncread( filename_help_fields,'v_3D',[1,1,ti],[Inf,Inf,1]);
  w_3D = ncread( filename_help_fields,'w_3D',[1,1,ti],[Inf,Inf,1]);

  trans.Hi = trans.A * Hi;

  trans.x_3D = zeros( trans.n, nz);
  trans.z_3D = zeros( trans.n, nz);
  trans.Ti   = zeros( trans.n, nz);
  trans.u_3D = zeros( trans.n, nz);
  trans.v_3D = zeros( trans.n, nz);
  trans.w_3D = zeros( trans.n, nz);

  for k = 1: nz
    trans.x_3D( :,k) = trans.x;
    trans.z_3D( :,k) = trans.Hi * (1 - zeta( k));
    trans.Ti(   :,k) = trans.A * Ti(   :,k);
    trans.u_3D( :,k) = trans.A * u_3D( :,k);
    trans.v_3D( :,k) = trans.A * v_3D( :,k);
    trans.w_3D( :,k) = trans.A * w_3D( :,k);
  end

  %% Update plot

  % Temperature
  [V,F,C] = xyc2vfc( trans.x_3D / 1e3, trans.z_3D, trans.Ti);
  set(H.Patch,'vertices',V,'faces',F,'facevertexcdata',C);

  % Velocities
  if ~isempty(Quiv); delete(Quiv); end
  hold on
  Quiv = quiver( trans.x_3D / 1e3, trans.z_3D, trans.u_3D / 1e3, trans.w_3D,'k');

  % Ice sheet surface
  set(H.IceSurf,'xdata',trans.x / 1e3,'ydata',trans.Hi);
  
  drawnow('update')
  pause(dt_pause)

end

function [V,F,C] = xyc2vfc( xx,yy,cc)

nx = size( xx,1);
ny = size( xx,2);

ij2n = zeros( nx,ny);

n = 0;
for i = 1: nx
  for j = 1: ny
    n = n+1;
    ij2n( i,j) = n;
  end
end

nV = nx*ny;
nF = (nx-1) * (ny-1);

V = zeros( nV,2);
F = zeros( nF,4);
C = zeros( nV,1);

for i = 1: nx
  for j = 1: ny
    n = ij2n( i,j);
    V( n,:) = [xx( i,j), yy( i,j)];
    C( n  ) = cc( i,j);
  end
end

fi = 0;
for i = 1: nx-1
  for j = 1: ny-1
    fi = fi+1;
    vi_sw = ij2n( i  ,j  );
    vi_se = ij2n( i+1,j  );
    vi_ne = ij2n( i+1,j+1);
    vi_nw = ij2n( i  ,j+1);
    F( fi,:) = [vi_sw, vi_se, vi_ne, vi_nw];
  end
end

end