clc
clear all
close all

experiments = {
  'historical'
%   'ctrl'
%   'NorESM1-M_RCP26-repeat'
%   'CCSM4_RCP85'
%   'HadGEM2-ES_RCP85'
%   'CESM2-WACCM_ssp585'
%   'UKESM1-0-LL_ssp585'
%   'UKESM1-0-LL_ssp585-repeat'
%   'NorESM1-M_RCP85-repeat'
%   'CCSM4_RCP85_shelfcollapse'
%   'HadGEM2-ES_RCP85_shelfcollapse'
%   'CESM2-WACCM_ssp585_shelfcollapse'
%   'UKESM1-0-LL_ssp585_shelfcollapse'
  };
experiment_codes = {
  'hist'
%   'ctrlAE'
%   'expAE01'
%   'expAE02'
%   'expAE03'
%   'expAE04'
%   'expAE05'
%   'expAE06'
%   'expAE07'
%   'expAE11'
%   'expAE12'
%   'expAE13'
%   'expAE14'
  };

model_versions = {
  'results_ANT_30km_ismip_v1'
  'results_ANT_30km_ismip_v2'
  'results_ANT_16km_ismip_v3'
  'results_ANT_16km_ismip_v4'
  };
model_versions_short = {
  'v1'
  'v2'
  'v3'
  'v4'
  };
model_names = {
  'UFEMISM1'
  'UFEMISM2'
  'UFEMISM3'
  'UFEMISM4'
  };

parts = {...
  '2014-2100'
  '2100-2300'
  };

for xi = 1: length( experiments)
  
  experiment      = experiments{      xi};
  experiment_code = experiment_codes{ xi};
  
  for mi = 1: length( model_versions)
    
    model_version       = model_versions{       mi};
    model_version_short = model_versions_short{ mi};
    model_name          = model_names{          mi};
    
    if strcmpi( experiment,'historical')
      parts_applied = {''};
    else
      parts_applied = parts;
    end
    
    for parti = 1: length( parts_applied)
      
      part = parts_applied{ parti};
    
      % Base UFEMISM results folder
      foldername_base = ['/Users/berends/Documents/Models/UFEMISM/' model_version '_' experiment '_' part];
      if strcmpi( foldername_base( end),'_')
        foldername_base = foldername_base( 1:end-1);
      end

      % Safety
      if ~exist( foldername_base,'dir')
        error('whaa!');
      end
      
      % ISMIP6 results folder
      foldername_ISMIP = [foldername_base '/AIS_IMAU_' model_name '_' experiment_code];
      
      % Set up Netcdf template
      filename_topg = [foldername_ISMIP '/topg_AIS_IMAU_' model_name '_' experiment_code '.nc'];
      f = ncinfo( filename_topg);
      
      f.Variables( 4).Name = 'base';
      f.Variables( 4).Attributes( 1).Value = 'base_altitude';
      
      f.Filename = [foldername_ISMIP '/base_AIS_IMAU_' model_name '_' experiment_code '.nc'];
      
      % Create NetCDF file
      if exist( f.Filename,'file')
        delete( f.Filename)
      end
      ncwriteschema( f.Filename, f);
      
      % Copy grid and time data
      x    = ncread( filename_topg,'x'   );
      y    = ncread( filename_topg,'y'   );
      time = ncread( filename_topg,'time');
      
      ncwrite( f.Filename,'x'   , x);
      ncwrite( f.Filename,'y'   , y);
      ncwrite( f.Filename,'time', time);
      
      %% Read, calculate, and write ice geometry data
      
      filename_restart = [foldername_base '/restart_ANT_00001.nc'];
      time_restart = ncread( filename_restart,'time');
      
      % Set up remapping stuff
      if     contains( filename_restart,'30km')
        M_map = read_CSR_from_NetCDF( '/Users/berends/Documents/Models/UFEMISM/M_map_mesh2grid_30km.nc');
      elseif contains( filename_restart,'16km')
        M_map = read_CSR_from_NetCDF( '/Users/berends/Documents/Models/UFEMISM/M_map_mesh2grid_16km.nc');
      else
        error('whaa!');
      end
      
      % Grid-to-vector translation tables
      nx = length( x);
      ny = length( y);
      n2ij = zeros( nx*ny,2);
      n = 0;
      for i = 1: nx
        if (mod(i,2) == 1)
          for j = 1: ny
            n = n+1;
            n2ij( n,:) = [i,j];
          end
        else
          for j = ny: -1: 1
            n = n+1;
            n2ij( n,:) = [i,j];
          end
        end
      end
      
      pprev = -Inf;
      
      for ti = 1: length( time)
        
        prog = floor( ti * 100 / length( time));
        if prog >= pprev + 10
          pprev = prog;
          disp(['Calculating ice base elevation for ' experiment_code '_' model_version '_' part ' - progress: ' num2str( prog) '%...'])
        end
        
        % Find timeframe to read
        tj = find( time_restart == time( ti)/360);
        
        % Safety
        if isempty( tj)
          if time( ti) == 2300*360 && time_restart( end) == 2299
            tj = length( time_restart);
          else
            error('whaa!')
          end
        elseif length( tj) == 2
          tj = 2;
        end
        
        % Read data
        Hi = ncread( filename_restart,'Hi',[1,tj],[Inf,1]);
        Hs = ncread( filename_restart,'Hs',[1,tj],[Inf,1]);
        
        % Calculate ice base elevation
        Hib = Hs - Hi;
        
        % Map to grid
        Hib_grid = map_mesh_to_grid( M_map, nx, ny, n2ij, Hib);
        
        % Write to NetCDF
        ncwrite( f.Filename,'base',Hib_grid,[1,1,ti]);
        
      end % for ti = 1: length( time)
    
    end
  end
end

function d_grid = map_mesh_to_grid( M_map, nx, ny, n2ij, d_mesh)

d_grid = zeros( nx, ny);

d_grid_vec = M_map * d_mesh;

% Reshape data from vector form to grid form
for n = 1: nx*ny
  i = n2ij( n,1);
  j = n2ij( n,2);
  d_grid( i,j) = d_grid_vec( n);
end

end