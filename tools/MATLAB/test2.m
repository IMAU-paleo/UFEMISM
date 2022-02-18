clc
clear all
close all

foldername = '../results_20211216_001';

A  = full(read_CSR_from_NetCDF( [foldername '/A.nc']))
B  = full(read_CSR_from_NetCDF( [foldername '/B.nc']))
AB = full(read_CSR_from_NetCDF( [foldername '/AB.nc']))