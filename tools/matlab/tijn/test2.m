clc
clear all
close all

foldernames = {
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_v1_historical'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_v1_CCSM4_RCP85_2014-2100'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_v1_CCSM4_RCP85_2100-2300'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_v2_historical'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_v2_CCSM4_RCP85_2014-2100'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_v2_CCSM4_RCP85_2100-2300'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_16km_ismip_v3_historical'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_16km_ismip_v3_CCSM4_RCP85_2014-2100'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_16km_ismip_v3_CCSM4_RCP85_2100-2300'
  };

foldernames = {
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_v1_historical'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_v1_CESM2-WACCM_ssp585_2014-2100'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_v1_CESM2-WACCM_ssp585_2100-2300'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_v2_historical'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_v2_CESM2-WACCM_ssp585_2014-2100'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_v2_CESM2-WACCM_ssp585_2100-2300'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_16km_ismip_v4_historical'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_16km_ismip_v4_CESM2-WACCM_ssp585_2014-2100'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_16km_ismip_v4_CESM2-WACCM_ssp585_2100-2300'
  };

ocean_area       = 3.611E14;
seawater_density = 1028;
ice_density      = 910;

% Old
for fi = 1: 3
  filename = [foldernames{fi} '/help_fields_ANT_00001.nc'];
  mesh = read_mesh_from_file( filename);
  time = ncread( filename,'time');
  VAF = zeros( size( time));
  for ti = 1: length( time)
    Hi = ncread( filename,'Hi',[1,ti],[Inf,1]);
    Hb = ncread( filename,'Hb',[1,ti],[Inf,1]);
    TAF = thickness_above_floatation( Hi, Hb, Hi*0);
    VAF( ti) = sum( TAF .* mesh.A) * ice_density / ocean_area / seawater_density;
  end
  if (fi == 1)
    VAF0 = VAF( 1);
  end
  VAF = VAF - VAF0;
  VAF( end) = VAF( end-1);
  plot( time,VAF,'color','r','linewidth',3)
  hold on
end

% New
for fi = 4: 6
  filename = [foldernames{fi} '/help_fields_ANT_00001.nc'];
  mesh = read_mesh_from_file( filename);
  time = ncread( filename,'time');
  VAF = zeros( size( time));
  for ti = 1: length( time)
    Hi = ncread( filename,'Hi',[1,ti],[Inf,1]);
    Hb = ncread( filename,'Hb',[1,ti],[Inf,1]);
    TAF = thickness_above_floatation( Hi, Hb, Hi*0);
    VAF( ti) = sum( TAF .* mesh.A) * ice_density / ocean_area / seawater_density;
  end
  if (fi == 4)
    VAF0 = VAF( 1);
  end
  VAF = VAF - VAF0;
  VAF( end) = VAF( end-1);
  plot( time,VAF,'color','g','linewidth',3)
  hold on
end

% New 16 km
for fi = 7: 9
  filename = [foldernames{fi} '/help_fields_ANT_00001.nc'];
  mesh = read_mesh_from_file( filename);
  time = ncread( filename,'time');
  VAF = zeros( size( time));
  for ti = 1: length( time)
    Hi = ncread( filename,'Hi',[1,ti],[Inf,1]);
    Hb = ncread( filename,'Hb',[1,ti],[Inf,1]);
    TAF = thickness_above_floatation( Hi, Hb, Hi*0);
    VAF( ti) = sum( TAF .* mesh.A) * ice_density / ocean_area / seawater_density;
  end
  if (fi == 7)
    VAF0 = VAF( 1);
  end
  VAF = VAF - VAF0;
  VAF( end) = VAF( end-1);
  plot( time,VAF,'color','b','linewidth',3)
  hold on
end

set(gca,'fontsize',24,'xgrid','on','ygrid','on');
xlabel('Year')
ylabel('\DeltaVAF (m.s.l.e.)')