clc
clear all
close all

foldernames = {
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_1950-2014_ocean1950-1979mean'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_1950-1995_oceanccsm4_rcp8.5'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_1995-2100_oceanccsm4_rcp8.5_from_histor'
  '/Users/berends/Documents/Models/UFEMISM/results_ANT_30km_ismip_1995-2100_oceanccsm4_rcp8.5'
  };

results = [];
for fi = 1: length( foldernames)
  
%   filename = [foldernames{ fi} '/restart_ANT_00001.nc'];
%   
%   time = ncread( filename,'time');
%   mesh = read_mesh_from_file( filename);
%   
%   results( fi).time = time;
%   results( fi).VAF  = time*0;
%   
%   for ti = 1: length( time)
%     Hi  = ncread( filename,'Hi',[1,ti],[Inf,1]);
%     Hb  = ncread( filename,'Hb',[1,ti],[Inf,1]);
%     seawater_density = 1014;
%     ice_density      = 910;
%     TAF = max( 0, Hi - max( 0, -Hb * seawater_density / ice_density));
%     results( fi).VAF( ti) = sum( TAF .* mesh.A) / 3.611E14;
%   end
  
  filename = [foldernames{ fi},'/scalar_output_ANT.nc'];
  results( fi).time = ncread( filename,'time');
  results( fi).VAF  = ncread( filename,'ice_volume_af');
  
  plot( results( fi).time, results( fi).VAF, 'linewidth',3); hold on
  
end

legend(...
  'Historical GCM ensemble mean',...
  'Historical CCSM4',...
  'Future CCSM4 from ensemble mean histor',...
  'Future CCSM4 from own histor',...
  'location','southwest');
set(gcf,'color','w');
set(gca,'fontsize',20,'xgrid','on','ygrid','on')
ylabel(gca,'Volume above flotation (m.s.l.e.)')
xlabel('Year')