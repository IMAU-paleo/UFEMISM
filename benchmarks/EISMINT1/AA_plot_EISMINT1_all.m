clc
clear all
close all

foldernames = {
  'EISMINT1_B_SIASSA';
  };

%% Read model output
for fi = 1: length( foldernames)
  
  filename_restart     = [foldernames{ fi} '/restart_ANT_00001.nc'    ];
  filename_help_fields = [foldernames{ fi} '/help_fields_ANT_00001.nc'];
  
  results( fi).time = ncread( filename_restart,'time');
  
  zeta = ncread( filename_restart,'zeta');
  nz = length( zeta);
  
  % Read mesh
  mesh = read_mesh_from_file( filename_restart);
  
  % Find triangle containing the ice divide
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
  Hi_a = ncread( filename_restart, 'Hi', [via,1],[1,Inf]);
  Hi_b = ncread( filename_restart, 'Hi', [vib,1],[1,Inf]);
  Hi_c = ncread( filename_restart, 'Hi', [vic,1],[1,Inf]);
  results( fi).Hi = wa * Hi_a + wb * Hi_b + wc * Hi_c;
  
  Ti_basal_a = ncread( filename_help_fields, 'Ti_basal', [via,1],[1,Inf]);
  Ti_basal_b = ncread( filename_help_fields, 'Ti_basal', [vib,1],[1,Inf]);
  Ti_basal_c = ncread( filename_help_fields, 'Ti_basal', [vic,1],[1,Inf]);
  results( fi).Ti_basal = wa * Ti_basal_a + wb * Ti_basal_b + wc * Ti_basal_c;
  
  Ti_pmp_a = permute( ncread( filename_help_fields, 'Ti_pmp', [via,nz,1],[1,1,Inf]), [1,3,2]);
  Ti_pmp_b = permute( ncread( filename_help_fields, 'Ti_pmp', [vib,nz,1],[1,1,Inf]), [1,3,2]);
  Ti_pmp_c = permute( ncread( filename_help_fields, 'Ti_pmp', [vic,nz,1],[1,1,Inf]), [1,3,2]);
  results( fi).Ti_pmp = wa * Ti_pmp_a + wb * Ti_pmp_b + wc * Ti_pmp_c;
  
  results( fi).Ti_basal_rel = results( fi).Ti_basal - results( fi).Ti_pmp;
  
end

%% Set up GUI