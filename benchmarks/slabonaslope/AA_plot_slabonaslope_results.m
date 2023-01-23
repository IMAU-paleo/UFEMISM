clc
clear all
close all

% Plot colors
colour_SIA  = [1.0,0.1,0.1];
colour_DIVA = [1.0,0.7,0.2];
colour_BPA  = [0.8,0.1,0.7];

% Model parameters
A           = 1e-16;        % Glen's flow law factor      [Pa-3 yr-1]
n           = 3;            % Glen's flow law exponent    [-]
Hi          = 2000;         % Ice thickness               [m]
Hb          = 0;            % Bedrock elevation           [m]
Hs          = Hi + Hb;      % Surface elevation           [m]
dh_dx       = -0.01;        % Surface slope               [m/m]
ice_density = 910;          % Ice density                 [kg m-3]
grav        = 9.81;         % Acceleration of gravity     [m s-2]

% Calculate the analytical solution on a high-resoution vertical grid for plotting
an.z = linspace( Hb, Hs, 1000);
[an.u, an.du_dz, an.eta] = calc_SIA_analytical( A, n, Hi, Hb, Hs, dh_dx, ice_density, grav, an.z);
figure; plot(an.z,an.u,'linewidth',2,'color','k');

%% Read UFEMISM output - SIA
foldername = 'slabonaslope_SIA';
filename = [foldername '/help_fields_ANT_00001.nc'];
mesh = read_mesh_from_file( filename);
A = calc_transect_matrix_a( mesh, 0, 0);
zeta   = ncread( filename,'zeta'); nz = length( zeta);
time   = ncread( filename,'time'); ti = length( time);
u_3D   = ncread( filename,'u_3D',[1,1,ti],[Inf,Inf,1]);
u_prof = zeros( nz,1);
for k = 1: nz
  u_prof( k) = A * u_3D( :,k);
end
line('xdata',Hi*(1-zeta),'ydata',u_prof,'linestyle','none','marker','o','markersize',8,'markerfacecolor',colour_SIA,'markeredgecolor',colour_SIA);

%% Read UFEMISM output - DIVA
foldername = 'slabonaslope_DIVA';
filename = [foldername '/help_fields_ANT_00001.nc'];
mesh = read_mesh_from_file( filename);
A = calc_transect_matrix_a( mesh, 0, 0);
zeta   = ncread( filename,'zeta'); nz = length( zeta);
time   = ncread( filename,'time'); ti = length( time);
u_3D   = ncread( filename,'u_3D',[1,1,ti],[Inf,Inf,1]);
u_prof = zeros( nz,1);
for k = 1: nz
  u_prof( k) = A * u_3D( :,k);
end
line('xdata',Hi*(1-zeta),'ydata',u_prof,'linestyle','none','marker','*','markersize',8,'markerfacecolor',colour_DIVA,'markeredgecolor',colour_DIVA);

%% Read UFEMISM output - BPA
foldername = 'slabonaslope_BPA';
filename = [foldername '/help_fields_ANT_00001.nc'];
mesh = read_mesh_from_file( filename);
A = calc_transect_matrix_a( mesh, 0, 0);
zeta   = ncread( filename,'zeta'); nz = length( zeta);
time   = ncread( filename,'time'); ti = length( time);
u_3D   = ncread( filename,'u_3D',[1,1,ti],[Inf,Inf,1]);
u_prof = zeros( nz,1);
for k = 1: nz
  u_prof( k) = A * u_3D( :,k);
end
line('xdata',Hi*(1-zeta),'ydata',u_prof,'linestyle','none','marker','^','markersize',8,'markerfacecolor',colour_BPA,'markeredgecolor',colour_BPA);

function [u,du_dz,eta] = calc_SIA_analytical( A, n, Hi, Hb, Hs, dh_dx, ice_density, grav, z)
% Calculate the analytical solution to the SIA for a constant flow factor A

% Calculate the integral term
% int_b_z (h-zeta)^n dzeta = 1/(n+1) [ H^(n+1) - (h-z)^(n+1)]
int_hminzn = 1/(n+1) * (Hi^(n+1) - (Hs-z).^(n+1));

% Calculate the velocity u
u = -2 * (ice_density * grav)^n * abs(dh_dx)^(n-1) * dh_dx * A * int_hminzn;

% Calculate the vertical shear strain rate du/dz
du_dz = -2.0 * (ice_density * grav)^n * abs(dh_dx)^(n-1) * dh_dx * A * (Hs - z).^n;

% Calculate the effective viscosity eta
eta = 0.5 * A^(-1.0/n) * (0.25 * du_dz.^2).^((1.0 - n)/(2.0 * n));

end