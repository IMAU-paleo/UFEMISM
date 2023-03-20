clc
clear all
close all

experiments = {
  'ctrl'
  'NorESM1-M_RCP26-repeat'
  'CCSM4_RCP85'
  'HadGEM2-ES_RCP85'
  'CESM2-WACCM_ssp585'
  'UKESM1-0-LL_ssp585'
  'UKESM1-0-LL_ssp585-repeat'
  'NorESM1-M_RCP85-repeat'
  'CCSM4_RCP85_shelfcollapse'
  'HadGEM2-ES_RCP85_shelfcollapse'
  'CESM2-WACCM_ssp585_shelfcollapse'
  'UKESM1-0-LL_ssp585_shelfcollapse'
  };
experiment_codes = {
  'ctrlAE'
  'expAE01'
  'expAE02'
  'expAE03'
  'expAE04'
  'expAE05'
  'expAE06'
  'expAE07'
  'expAE11'
  'expAE12'
  'expAE13'
  'expAE14'
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


for xi = 1: length( experiments)
  
  experiment      = experiments{      xi};
  experiment_code = experiment_codes{ xi};
  
  for mi = 1: length( model_versions)
    
    model_version       = model_versions{       mi};
    model_version_short = model_versions_short{ mi};
    model_name          = model_names{          mi};
    
    foldername_src1 = ['/Users/berends/Documents/Models/UFEMISM/' model_version '_' experiment '_2014-2100/AIS_IMAU_' model_name '_' experiment_code];
    foldername_src2 = ['/Users/berends/Documents/Models/UFEMISM/' model_version '_' experiment '_2100-2300/AIS_IMAU_' model_name '_' experiment_code];
    
    foldername_dst  = ['/Users/berends/Documents/Models/UFEMISM/AIS_IMAU_' model_name '/AIS_IMAU_' model_name '_' experiment_code];
    
    merge_ISMIP_folder( foldername_src1, foldername_src2, foldername_dst);
    
  end % for mi = 1: length( model_versions)
end % for xi = 1: length( experiments)

function merge_ISMIP_folder( foldername_src1, foldername_src2, foldername_dst)

% Safety
if ~exist( foldername_src1,'dir')
  error('whaa!')
end
if ~exist( foldername_src2,'dir')
  error('whaa!')
end
if exist( foldername_dst,'dir')
  error('whaa!')
end

% Make output folder
mkdir( foldername_dst)

%% Safety
% Check if each file in dir1 also exists in dir2
dir_src1 = dir( foldername_src1);
dir_src2 = dir( foldername_src2);

for i = 1: length( dir_src1)
  if strcmpi( dir_src1( i).name,'.DS_store'); continue; end
  exists_in_other = false;
  for j = 1: length( dir_src2)
    if strcmpi( dir_src1( i).name, dir_src2( j).name)
      exists_in_other = true;
    end
  end
  if ~exists_in_other
    error('whaa!')
  end
end
for j = 1: length( dir_src2)
  if strcmpi( dir_src2( j).name,'.DS_store'); continue; end
  exists_in_other = false;
  for i = 1: length( dir_src1)
    if strcmpi( dir_src1( i).name, dir_src2( j).name)
      exists_in_other = true;
    end
  end
  if ~exists_in_other
    error('whaa!')
  end
end

disp(['Creating merged ISMIP folder "' foldername_dst '"...'])

%% Merge files
for i = 1: length( dir_src1)
  
  if ~contains( dir_src1( i).name,'.nc'); continue; end
  
  % Skip GHF as it has no time dimension
  if contains( dir_src1( i).name,'hfgeoubed'); continue; end
  
  filename_src1 = [foldername_src1 '/' dir_src1( i).name];
  filename_src2 = [foldername_src2 '/' dir_src1( i).name];
  filename_dst  = [foldername_dst  '/' dir_src1( i).name];
  
  disp(['Creating merged file "' filename_dst '"...'])
  
  % Read source times
  time_src1 = ncread( filename_src1,'time');
  time_src2 = ncread( filename_src2,'time');
  
  % Safety
  if time_src1( end) ~= time_src2( 1)
    error('whaa!')
  end
  
  % Remove double timeframes if needed
  ti_src1 = (1: length( time_src1))';
  while time_src1( 1) == time_src1( 2)
    ti_src1(   1) = [];
    time_src1( 1) = [];
  end
  
  ti_src2 = (1: length( time_src2))';
  while time_src2( 1) == time_src2( 2)
    ti_src2(   1) = [];
    time_src2( 1) = [];
  end
  
  % Merge time variable
  time_dst = [time_src1( 1:end-1); time_src2];
  ti_src1p = [ti_src1(   1:end-1); ti_src2*0];
  ti_src2p = [ti_src1(   1:end-1)*0; ti_src2];
  
  % Netcdf template
  f = ncinfo( filename_src1);
  
  di_time = 0;
  for di = 1: length( f.Dimensions)
    if strcmpi( f.Dimensions( di).Name,'time')
      di_time = di;
      f.Dimensions( di).Length = length( time_dst);
    end
  end
  % Safety
  if di_time == 0
    error('whaa!')
  end
  
  var_type = 'scalar';
  var_name = '';
  
  for vi = 1: length( f.Variables)
    if     strcmpi( f.Variables( vi).Name,'x')
      nx = f.Variables( vi).Size;
      var_type = 'field';
    elseif strcmpi( f.Variables( vi).Name,'y')
      ny = f.Variables( vi).Size;
      var_type = 'field';
    elseif strcmpi( f.Variables( vi).Name,'time')
      f.Variables( vi).Dimensions = f.Dimensions( di_time);
      f.Variables( vi).Size = length( time_dst);
    else
      var_name = f.Variables( vi).Name;
      found_time = false;
      for di = 1: length( f.Variables( vi).Dimensions)
        if strcmpi( f.Variables( vi).Dimensions( di).Name,'time')
          found_time = true;
          f.Variables( vi).Dimensions( di).Length = length( time_dst);
          f.Variables( vi).Size( di) = length( time_dst);
        end
      end
      % Safety
      if ~found_time
        error('whaa!')
      end
    end
  end
  
  % Create merged file
  ncwriteschema( filename_dst,f);
  
  % Copy data
  if strcmpi( var_type,'scalar')
    
    % Write time
    ncwrite( filename_dst, 'time', time_dst);
    
    % Memory for merged data
    d_dst = zeros( length( time_dst),1);
    
    % Read data from both source files
    for tii = 1: length( time_dst)
      ti_src1 = ti_src1p( tii);
      ti_src2 = ti_src2p( tii);
      if ti_src1 > 0
        d_dst( tii) = ncread( filename_src1, var_name, ti_src1, 1);
      else
        d_dst( tii) = ncread( filename_src2, var_name, ti_src2, 1);
      end
    end
    
    % Write to destination file
    ncwrite( filename_dst, var_name, d_dst);
    
  else
    
    % Write x,y,time
    ncwrite( filename_dst, 'x'   , ncread( filename_src1, 'x'));
    ncwrite( filename_dst, 'y'   , ncread( filename_src1, 'y'));
    ncwrite( filename_dst, 'time', time_dst);
    
    % Memory for merged data
    d_dst = zeros( nx, ny, length( time_dst),1);
    
    % Read data from both source files
    for tii = 1: length( time_dst)
      ti_src1 = ti_src1p( tii);
      ti_src2 = ti_src2p( tii);
      if ti_src1 > 0
        d_dst( :,:,tii) = ncread( filename_src1, var_name, [1,1,ti_src1], [Inf,Inf,1]);
      else
        d_dst( :,:,tii) = ncread( filename_src2, var_name, [1,1,ti_src2], [Inf,Inf,1]);
      end
    end
    
    % Write to destination file
    ncwrite( filename_dst, var_name, d_dst);
    
  end
  
end

end