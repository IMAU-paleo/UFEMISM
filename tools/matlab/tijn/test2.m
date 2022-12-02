clc
clear all
close all

foldername = '../../../benchmarks/EISMINT1/EISMINT1_A';
filename   = [foldername '/debug_ANT.nc'];
mesh = read_mesh_from_file( filename);

mask_ice_a           = ncread( filename, 'int_2D_a_01');
Hi_a                 = ncread( filename, 'dp_2D_a_01');
u_3D_a               = ncread( filename, 'dp_3D_a_01');
v_3D_a               = ncread( filename, 'dp_3D_a_02');
w_3D_a               = ncread( filename, 'dp_3D_a_03');
frictional_heating_a = ncread( filename, 'dp_2D_a_02');
internal_heating_a   = ncread( filename, 'dp_3D_a_04');
T_surf_annual        = ncread( filename, 'dp_2D_a_03');
Q_base_grnd          = ncread( filename, 'dp_2D_a_04');
T_base_float         = ncread( filename, 'dp_2D_a_05');
Ti_a                 = ncread( filename, 'dp_3D_a_05');
Ti_tplusdt           = ncread( filename, 'dp_3D_a_06');
is_unstable          = ncread( filename, 'int_2D_a_02');

plot_mesh_data( mesh, is_unstable, 'k');