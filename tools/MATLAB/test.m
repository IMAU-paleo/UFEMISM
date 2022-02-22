clc
clear all
close all

foldername = '../results_20211216_001';

mesh = ReadMeshFromFile( [foldername '/restart_ANT_00001.nc']);

f_c  = ncread( [foldername '/debug_ANT.nc'],'dp_2D_c_01');

PlotMeshData( mesh, f_c, 'k');