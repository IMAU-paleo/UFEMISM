clc
clear all
close all

%% In/output folder names
foldername_src1 = '/Users/berends/Documents/Models/UFEMISM/results_old/results_ANT_30km_ismip_v1_CCSM4_RCP85_2014-2100/AIS_IMAU_UFEMISM_histor';
foldername_src2 = '/Users/berends/Documents/Models/UFEMISM/results_old/results_ANT_30km_ismip_v1_CCSM4_RCP85_2100-2300/AIS_IMAU_UFEMISM_histor';

foldername_dst  = '/Users/berends/Documents/Models/UFEMISM/AIS_IMAU_UFEMISM_histor';

%% Clean output folder
if exist(  foldername_dst,'dir')
  delete( [foldername_dst '/*'])
  rmdir(   foldername_dst)
end
mkdir( foldername_dst)

%% Safety
dir_src1 = dir( foldername_src1);
dir_src2 = dir( foldername_src2);

for i = 1: length( dir_src1)
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

%% Merge files
for i = 1: length( dir_src1)
  
  if ~contains( dir_src1( i).name,'.nc'); continue; end
  
  % Skip GHF as it has no time dimension
  if contains( dir_src1( i).name,'hfgeoubed'); continue; end
  
  filename_src1 = [foldername_src1 '/' dir_src1( i).name];
  filename_src2 = [foldername_src2 '/' dir_src1( i).name];
  filename_dst  = [foldername_dst  '/' dir_src1( i).name];
  
  time_src1 = ncread( filename_src1,'time');
  time_src2 = ncread( filename_src2,'time');
  
  % Safety
  if time_src1( end) ~= time_src2( 1)
    error('whaa!')
  end
  
  time_dst = [time_src1( 1:end-1); time_src2];
  
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
  if di_time == 0; error('whaa!'); end
  
  var_type = 'scalar';
  var_name = '';
  
  for vi = 1: length( f.Variables)
    if     strcmpi( f.Variables( vi).Name,'x')
      var_type = 'field';
    elseif strcmpi( f.Variables( vi).Name,'y')
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
      if ~found_time; error('whaa!'); end
    end
  end
  
  % Create merged file
  ncwriteschema( filename_dst,f);
  
  % Copy data
  if strcmpi( var_type,'scalar')
    d_src1 = ncread( filename_src1, var_name, 1, length( time_src1)-1);
    d_src2 = ncread( filename_src2, var_name);
    d_dst  = [d_src1; d_src2];
    ncwrite( filename_dst, var_name, d_dst);
  else
    d_src1 = ncread( filename_src1, var_name, [1,1,1], [Inf,Inf,length( time_src1)-1]);
    d_src2 = ncread( filename_src2, var_name);
    d_dst  = cat( 3, d_src1, d_src2);
    ncwrite( filename_dst, var_name, d_dst);
  end
  
end