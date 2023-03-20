clc
clear all
close all

%% Initialise basic GUI

wa          = 700;
ha          = 500;
margins_hor = [100,25];
margins_ver = [25,100];
ylim        = [0,8];

H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

set(    H.Ax{1,1},'xlim',[0,500],'ylim',ylim);
xlabel( H.Ax{1,1},'Year')
ylabel( H.Ax{1,1},'Sea level contribution (m)')

%% Experiments to plot

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

% Parameters
seawater_density = 1028;
ice_density      = 910;
ocean_area       = 3.611e14;

  
for mi = 1: length( model_versions)

  model_version       = model_versions{       mi};
  model_version_short = model_versions_short{ mi};

  disp(['Reading results for abumip, model version ' model_version_short '...'])

  %% Part 1

  % UFEMISM results folder
  foldername = ['/Users/berends/Documents/Models/UFEMISM/' model_version '_abumip'];

  filename = [foldername '/restart_ANT_00001.nc'];
  
  mesh = read_mesh_from_file( filename);

  % Read data
  time = ncread( filename,'time');
  VAF_MSLE = time*0;
  for ti = 1: length( time)
    Hi = ncread( filename,'Hi',[1,ti],[Inf,1]);
    Hb = ncread( filename,'Hb',[1,ti],[Inf,1]);
    TAF = thickness_above_floatation( Hi, Hb, Hb*0);
    VAF = sum( max( TAF,0) .* mesh.A);
    VAF_MSLE( ti) = VAF * ice_density / (ocean_area * seawater_density);
  end
  
  VAF_MSLE_rel = VAF_MSLE - VAF_MSLE( 1);
  
  results.(model_version_short).time         = time(         2:end);
  results.(model_version_short).VAF_MSLE_rel = VAF_MSLE_rel( 2:end);

end % for mi = 1: length( model_versions)
  
%% Create plumes
  
results.time         = results.v1.time;
results.VAF_MSLE_min = results.v1.time * 0 + Inf;
results.VAF_MSLE_max = results.v1.time * 0 - Inf;
results.VAF_MSLE_av  = results.v1.time * 0;

for ti = 1: length( results.time)

  n = 0;

  for mi = 1: length( model_versions)

    model_version       = model_versions{       mi};
    model_version_short = model_versions_short{ mi};

    tj = find( results.(model_version_short).time == results.time( ti));

    if tj > 0
      n = n+1;
      VAF_MSLE = results.(model_version_short).VAF_MSLE_rel( tj);
      results.VAF_MSLE_min( ti) = min( results.VAF_MSLE_min( ti),  VAF_MSLE);
      results.VAF_MSLE_max( ti) = max( results.VAF_MSLE_max( ti),  VAF_MSLE);
      results.VAF_MSLE_av(  ti) =      results.VAF_MSLE_av(  ti) + VAF_MSLE;
    end

  end % for mi = 1: length( model_versions)

  results.VAF_MSLE_av( ti) = results.VAF_MSLE_av( ti) / n;

end % for ti = 1: length( results.(experiment_code).time)

%% Plot

color = 'r';

xdata = [results.time        ; flipud( results.time        )];
ydata = [results.VAF_MSLE_min; flipud( results.VAF_MSLE_max)];
patch('parent',H.Ax{1,1},'xdata',xdata,'ydata',-ydata,'facecolor',color,'edgecolor','none','facealpha',0.5);
line( 'parent',H.Ax{1,1},'xdata',results.time,'ydata',-results.VAF_MSLE_min,'color',color,'linewidth',2);
line( 'parent',H.Ax{1,1},'xdata',results.time,'ydata',-results.VAF_MSLE_max,'color',color,'linewidth',2);
line( 'parent',H.Ax{1,1},'xdata',results.time,'ydata',-results.VAF_MSLE_av ,'color',color,'linewidth',3);
  