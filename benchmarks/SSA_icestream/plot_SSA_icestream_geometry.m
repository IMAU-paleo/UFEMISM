clc
clear all
close all

H.Fig = figure('color','w');
H.Ax  = axes('color','w');

% Model parameters
C.A    = 1e-18;       % Glen's flow law factor [Pa^-3 yr^-1]
C.n    = 3;           % Glen's flow law exponent
C.H    = 2000;        % Ice thickness (slab on a sloe)
C.dhdx = 0.0003;      % Slope
C.m    = 1;           % Ice-stream friction exponent (Schoof)
C.L    = 150e3;       % Ice-stream half-width [m]
C.rho  = 910;         % Ice density [kg m^-3]
C.g    = 9.81;        % Acceleration of gravity [m s^-1]

an.y = linspace( -400e3, 400e3, 16);
an.u = SSA_Schoof2006_analytical_solution( C, an.y);

% Plot analytical solution
x = zeros( size( an.y)) + 0.25;
y = linspace( 0,1,length( an.y));
z = zeros( size( an.y)) + 2.75;
u = an.u;
v = zeros( size( an.y));
w = -an.u;
quiver3(H.Ax,x,y,z,u,v,w,3,'color','k','linewidth',2); hold on
quiver3(H.Ax,x,y,z-2,u,v,w,3,'color','k','linewidth',2)

% Plot ice slab

V = [
  0,0,0; 
  1,0,0;
  1,1,0;
  0,1,0;
  0,0,1;
  0,1,1;
  0,0,3;
  1,0,2;
  1,1,2;
  0,1,3];
F_bed = [
  1,2,3,4;
  2,3,6,5;
  1,4,6,5];
F_ice = [
  5,2,8,7;
  2,3,9,8;
  6,3,9,10;
  5,6,10,7;
  7,8,9,10];

set(H.Ax,'cameraposition',[6.8017   -4.7131   10.0443]);
H.Ax.XAxis.Visible = 'off';
H.Ax.YAxis.Visible = 'off';
H.Ax.ZAxis.Visible = 'off';
H.Patch_ice = patch('parent',H.Ax,'vertices',V,'faces',F_ice,'facecolor','none',...
  'edgecolor',[0.3,0.7,1.0],'linewidth',3);
H.Patch_bed = patch('parent',H.Ax,'vertices',V,'faces',F_bed,'facecolor','none',...
  'edgecolor',[0.8,0.6,0.1],'linewidth',3);

function u = SSA_Schoof2006_analytical_solution( C, y)

% Inverse flow factor (ce stiffness?)
B = C.A^(-1/C.n);

% Calculate the gravitational driving stress f
f = C.rho * C.g * C.H * C.dhdx;

% Calculate u
ua = -2 * f.^3 * C.L.^4 / (B.^3 * C.H.^3);
ub = ( 1 / 4                       ) * (   (y/C.L).^       4  - (C.m+1).^(   4/C.m) );
uc = (-3 / ((C.m+1)    * (  C.m+4))) * (abs(y/C.L).^(  C.m+4) - (C.m+1).^(1+(4/C.m)));
ud = ( 3 / ((C.m+1).^2 * (2*C.m+4))) * (abs(y/C.L).^(2*C.m+4) - (C.m+1).^(2+(4/C.m)));
ue = (-1 / ((C.m+1).^3 * (3*C.m+4))) * (abs(y/C.L).^(3*C.m+4) - (C.m+1).^(3+(4/C.m)));
u = ua * (ub + uc + ud + ue);

% Outside the ice-stream, velocity is zero
W = C.L * (C.m+1)^(1/C.m);
u( abs(y) > W) = 0;

end