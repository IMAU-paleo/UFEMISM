clc
clear all
close all

experiments = {
%   'ctrl'
%   'NorESM1-M_RCP26-repeat'
%   'CCSM4_RCP85'
%   'HadGEM2-ES_RCP85'
%   'CESM2-WACCM_ssp585'
%   'UKESM1-0-LL_ssp585'
  'UKESM1-0-LL_ssp585-repeat'
%   'NorESM1-M_RCP85-repeat'
%   'CCSM4_RCP85_shelfcollapse'
%   'HadGEM2-ES_RCP85_shelfcollapse'
%   'CESM2-WACCM_ssp585_shelfcollapse'
%   'UKESM1-0-LL_ssp585_shelfcollapse'
  };
experiment_codes = {
%   'ctrlAE'
%   'expAE01'
%   'expAE02'
%   'expAE03'
%   'expAE04'
%   'expAE05'
  'expAE06'
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

parts = {
  '2100-2300'
  };


for xi = 1: length( experiments)
  
  experiment      = experiments{      xi};
  experiment_code = experiment_codes{ xi};
  
  for mi = 1: length( model_versions)
    
    model_version       = model_versions{       mi};
    model_version_short = model_versions_short{ mi};
    model_name          = model_names{          mi};
    
    for partsi = 1: length( parts)
      
      part = parts{ partsi};
    
      % UFEMISM results folder
      foldername_base = ['/Users/berends/Documents/Models/UFEMISM/' model_version '_' experiment '_' part];
      
      % ISMIP6 results folder
      foldername_ISMIP = [foldername_base '/AIS_IMAU_' model_name '_' experiment_code];

      % Find all netcdf files in this folder
      netcdf_files = [];
      netcdf_files = list_all_netcdf_files_below( netcdf_files, foldername_ISMIP);

      % Add the last year to all of them
      for j = 1: length( netcdf_files)

        filename = netcdf_files{ j};

        % Check if this file has a time dimension
        has_time = false;
        f = ncinfo( filename);
        for di = 1: length( f.Dimensions)
          if strcmpi( f.Dimensions( di).Name,'time')
            has_time = true;
            break
          end
        end
        if ~has_time; continue; end

        % Add the last year
        add_last_year( filename,f);
      end
    
    end
  end
end

function netcdf_files = list_all_netcdf_files_below( netcdf_files, current_dir)
% List all netcdf files int he current folder, and in all subdirectories
% (recursively)

henk = dir( current_dir);

for i = 1: length( henk)
  if ~henk( i).isdir && contains( henk( i).name, '.nc')
    if contains( henk( i).name,'GCM_ensemble_mean'); continue; end
    netcdf_files{ end+1} = [current_dir '/' henk( i).name];
  elseif henk( i).isdir && ~(strcmpi( henk( i).name, '.') || strcmpi( henk( i).name, '..'))
    netcdf_files = list_all_netcdf_files_below( netcdf_files, [current_dir '/' henk( i).name]);
  end
end

end
function add_last_year( filename,f)
% Add the year 2300 to this file by copying the year 2299

disp(['Inspecting file "' filename '"...'])

% Safety
if  contains( filename,'resource_tracking.nc') || ...
    contains( filename,'restart') || ...
    contains( filename,'help_fields') || ...
    contains( filename,'scalar_output')
  return
end

% % Safety!
% filename_copy = strrep( filename,'.nc','_old.nc');
% copyfile( filename, filename_copy);

% Read time
time = ncread( filename,'time');

if length( time)<=1
  error('whaa!')
end

% If this file contains a 2300 frame, no need to do anything
if (time( end) == 2300 || time( end) == 828000)
  return
end

% Safety
if ~(time( end) == 2299 || time( end) == 827640)
  error('whaa!')
end

disp(['Appending year 2300 to file "' filename '"...'])

% Append 2300
ti_2299 = length( time);
time_2299 = time( ti_2299);
dt = time( ti_2299) - time( ti_2299-1);
time_2300 = time_2299 + dt;
time = [time; time_2300];

% Write to file
ncwrite( filename,'time',time);

% Copy year 2299 to year 2300 for all variables
for vi = 1: length( f.Variables)
  
  var = f.Variables( vi);
  
  % We already covered the time variable
  if strcmpi( var.Name,'time'); continue; end
  
  % Check if this variable has a time dimension; if not, skip it
  has_time = false;
  for di = 1: length( var.Dimensions)
    if strcmpi( var.Dimensions( di).Name,'time')
      has_time = true;
    end
  end
  if ~has_time; continue; end
  
  % Check if this is a scalar or a field variable
  if     length( var.Size) == 1
    var_type = 'scalar';
  elseif length( var.Size) == 2
    var_type = 'field_mesh';
  elseif length( var.Size) == 3
    var_type = 'field';
  else
    error('whaa!')
  end
  
  if strcmpi( var_type,'scalar')
    % Scalar variable
  
    % Read the year 2299
    d_2299 = ncread( filename, var.Name, ti_2299,1);
    
    % Append the year 2300
    ncwrite( filename, var.Name, d_2299, ti_2299+1);
    
  elseif strcmpi( var_type,'field_mesh')
    % Mesh field variable
  
    % Read the year 2299
    d_2299 = ncread( filename, var.Name, [1,ti_2299],[Inf,1]);
    
    % Append the year 2300
    ncwrite( filename, var.Name, d_2299, [1,ti_2299+1]);
    
  elseif strcmpi( var_type,'field')
    % Field variable
  
    % Read the year 2299
    d_2299 = ncread( filename, var.Name, [1,1,ti_2299],[Inf,Inf,1]);
    
    % Append the year 2300
    ncwrite( filename, var.Name, d_2299, [1,1,ti_2299+1]);
    
  else
    error('whaa!')
  end
end

end