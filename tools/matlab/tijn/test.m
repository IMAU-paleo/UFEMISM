clc
clear all
close all

filename1 = '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_1950-1995_oceanccsm4_rcp8.5/restart_ANT_00001.nc';
filename2 = '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_1995-2100_oceanccsm4_rcp8.5/restart_ANT_00001.nc';

mesh = read_mesh_from_file( filename1);

time1 = ncread( filename1,'time');
Hi1   = ncread( filename1,'Hi',[1,length( time1)],[Inf,1]);
Hb1   = ncread( filename1,'Hi',[1,length( time1)],[Inf,1]);

Hi2   = ncread( filename2,'Hi',[1,1],[Inf,1]);
Hb2   = ncread( filename2,'Hi',[1,1],[Inf,1]);

plot_mesh_data( mesh, Hi2 - Hi1,'none');
plot_mesh_data( mesh, Hb2 - Hb1,'none');