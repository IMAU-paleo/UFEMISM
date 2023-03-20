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
      
      % Correct name of ISMIP results folder
      foldername_ISMIP_correct = ['AIS_IMAU_' model_name '_' experiment_code];
      
      % Actual name of ISMIP results folder
      foldername_ISMIP_actual = [];
      henk = dir( foldername_base);
      for i = 1: length( henk)
        if henk( i).isdir && contains( henk( i).name,'AIS_IMAU')
          foldername_ISMIP_actual = henk( i).name;
          break
        end
      end
      % Safety
      if isempty( foldername_ISMIP_actual)
        error('whaa!')
      end
      
      % Rename
      if ~strcmpi( foldername_ISMIP_actual, foldername_ISMIP_correct)
        movefile( [foldername_base '/' foldername_ISMIP_actual], [foldername_base '/' foldername_ISMIP_correct])
      end
      
      foldername_ISMIP = [foldername_base '/' foldername_ISMIP_correct];
      
      % Individual files
      henk = dir( foldername_ISMIP);
      for i = 1: length( henk)
        if contains( henk( i).name,'AIS_IMAU')
          
          filename_actual = henk( i).name;
          j = strfind( filename_actual,'AIS_IMAU');
          % Safety
          if j==0
            error('whaa!')
          end
          varcode = filename_actual( 1:j-2);
          
          filename_correct = [varcode '_AIS_IMAU_' model_name '_' experiment_code '.nc'];
          
          % Rename
          if ~strcmpi( filename_actual, filename_correct)
            disp( [ 'Renaming "' filename_actual '" to "' filename_correct '"...'])
            movefile( [foldername_ISMIP '/' filename_actual], [foldername_ISMIP '/' filename_correct])
          end
          
        end
      end
      
      
    end % for parti = 1: 2
  end % for mi = 1: length( model_versions)
end % for xi = 1: length( experiments)