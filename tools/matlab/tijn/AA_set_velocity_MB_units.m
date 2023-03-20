clc
clear all
close all

experiments = {
  'historical'
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
  'hist'
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

parts = {...
  '2014-2100'
  '2100-2300'
  };

sec_per_year = 365 * 24 * 3600;
ice_density  = 910;

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
      
      disp(['Checking run ' foldername_base '...'])
      
      % ISMIP6 results folder
      foldername_ISMIP = [foldername_base '/AIS_IMAU_' model_name '_' experiment_code];
      
      henk = dir( foldername_ISMIP);
      for i = 1: length( henk)
        
        filename = henk( i).name;
        j = strfind( filename,'_AIS_IMAU');
        if isempty( j); continue; end
        var_name = filename( 1:j-1);
        
        if  strcmpi( var_name,'xvelbase') || ...
            strcmpi( var_name,'xvelmean') || ...
            strcmpi( var_name,'xvelsurf') || ...
            strcmpi( var_name,'yvelbase') || ...
            strcmpi( var_name,'yvelmean') || ...
            strcmpi( var_name,'yvelsurf') || ...
            strcmpi( var_name,'zvelbase') || ...
            strcmpi( var_name,'zvelsurf')
          % Velocities should be expressded in m/s, not in m/yr
          
          d = ncread( [foldername_ISMIP '/' filename], var_name);
          umax = max( abs( d(:)));
          
          if umax > 0.1
            % Yeah this is m/yr
            disp(['Fixing velocity ' var_name ' units in run ' foldername_base '...'])
            d = d / sec_per_year;
            ncwrite( [foldername_ISMIP '/' filename], var_name, d);
          end
          
        elseif strcmpi( var_name,'acabf') || ...
            strcmpi( var_name,'libmassbfgr') || ...
            strcmpi( var_name,'libmassbffl')
          % Mass balance components should be expressed in kg/m2/s, not in m.i.e./yr
          
          d = ncread( [foldername_ISMIP '/' filename], var_name);
          MB_max = max( abs( d(:)));
          
          if MB_max > 0.1
            % Yeah this is m.i.e./yr
            disp(['Fixing MB component ' var_name ' units in run ' foldername_base '...'])
            d = d * ice_density / sec_per_year;
            ncwrite( [foldername_ISMIP '/' filename], var_name, d);
          end
          
        else
          % No need to do anything to this variable
        end
        
      end
    
    end
  end
end