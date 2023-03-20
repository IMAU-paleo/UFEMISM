clc
clear all
close all

%% Initialise basic GUI

wa          = [700,50];
ha          = 500;
margins_hor = [100,5,100];
margins_ver = [25,100];
ylim        = [-2.0,1.0];

H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

set( H.Ax{1,2},'xgrid','off','ygrid','off','yaxislocation','right');
H.Ax{1,2}.XAxis.Visible = 'off';

set(    H.Ax{1,1},'xlim',[1950,2300],'ylim',ylim);
set(    H.Ax{1,2},'xlim',[0,1],'ylim',ylim);
xlabel( H.Ax{1,1},'Year')
ylabel( H.Ax{1,1},'Sea level contribution (m)')
ylabel( H.Ax{1,2},'Sea level contribution in 2100 (m)')

%% Experiments to plot

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

% Parameters
ice_density      = 910;
seawater_density = 1028;
ocean_area       = 3.611e14;

%% Read results for historical period
for mi = 1: length( model_versions)
  
  model_version       = model_versions{       mi};
  model_version_short = model_versions_short{ mi};
  
  disp(['Reading results for historical period, model version ' model_version_short '...'])
  
  % UFEMISM results folder
  foldername_base = ['/Users/berends/Documents/Models/UFEMISM/' model_version '_historical'];
  
  % Find the ISMIP6 output folder
  foldername_ISMIP = '';
  henk = dir( foldername_base);
  for i = 1: length( henk)
    if henk( i).isdir && contains( henk( i).name,'AIS_IMAU_UFEMISM')
      foldername_ISMIP = [foldername_base '/' henk( i).name];
      break
    end
  end
  % Safety
  if isempty( foldername_ISMIP)
    error(['Couldnt find ISMIP output folder in UFEMISM results folder "' foldername_base '"!'])
  end
  
  % Find the ISMIP6 output file for total grounded ice mass
  filename_ISMIP = '';
  henk = dir( foldername_ISMIP);
  for i = 1: length( henk)
    if contains( henk( i).name,'limnsw')
      filename_ISMIP = [foldername_ISMIP '/' henk( i).name];
    end
  end
  % Safety
  if isempty( filename_ISMIP)
    error(['Couldnt find limnsw file in ISMIP results folder "' foldername_ISMIP '"!'])
  end

  % Read data
  time     = ncread( filename_ISMIP,'time') / 360; % 360-day calendar
  VAF_MSLE = ncread( filename_ISMIP,'limnsw') / (seawater_density * ocean_area); % expressed in m.s.l.e.

  % DENK DROM - should fix this in the files, first time frame is written twice
  time     = time(     2:end);
  VAF_MSLE = VAF_MSLE( 2:end);
  
  % Save to results structure
  results.historical.(model_version_short).time         = time;
  results.historical.(model_version_short).VAF_MSLE     = VAF_MSLE;
  results.historical.(model_version_short).VAF_MSLE_rel = VAF_MSLE - VAF_MSLE( end);
  
end % for mi = 1: length( model_versions)

%% Read results for different experiments
for xi = 1: length( experiments)
  
  experiment      = experiments{      xi};
  experiment_code = experiment_codes{ xi};
  
  for mi = 1: length( model_versions)

    model_version       = model_versions{       mi};
    model_version_short = model_versions_short{ mi};

    disp(['Reading results for experiment ' experiment ', model version ' model_version_short '...'])
    
    %% Part 1

    % UFEMISM results folder
    foldername_base = ['/Users/berends/Documents/Models/UFEMISM/' model_version '_' experiment '_2014-2100'];

    % Find the ISMIP6 output folder
    foldername_ISMIP = '';
    henk = dir( foldername_base);
    for i = 1: length( henk)
      if henk( i).isdir && contains( henk( i).name,'AIS_IMAU_UFEMISM')
        foldername_ISMIP = [foldername_base '/' henk( i).name];
        break
      end
    end
    % Safety
    if isempty( foldername_ISMIP)
      error(['Couldnt find ISMIP output folder in UFEMISM results folder "' foldername_base '"!'])
    end

    % Find the ISMIP6 output file for total grounded ice mass
    filename_ISMIP = '';
    henk = dir( foldername_ISMIP);
    for i = 1: length( henk)
      if contains( henk( i).name,'limnsw')
        filename_ISMIP = [foldername_ISMIP '/' henk( i).name];
      end
    end
    % Safety
    if isempty( filename_ISMIP)
      error(['Couldnt find limnsw file in ISMIP results folder "' foldername_ISMIP '"!'])
    end

    % Read data
    time1     = ncread( filename_ISMIP,'time') / 360; % 360-day calendar
    VAF_MSLE1 = ncread( filename_ISMIP,'limnsw') / (seawater_density * ocean_area); % expressed in m.s.l.e.

    % DENK DROM - should fix this in the files, first time frame is written twice
    time1     = time1(     2:end);
    VAF_MSLE1 = VAF_MSLE1( 2:end);
    
    %% Part 2

    % UFEMISM results folder
    foldername_base = ['/Users/berends/Documents/Models/UFEMISM/' model_version '_' experiment '_2100-2300'];

    % Find the ISMIP6 output folder
    foldername_ISMIP = '';
    henk = dir( foldername_base);
    for i = 1: length( henk)
      if henk( i).isdir && contains( henk( i).name,'AIS_IMAU_UFEMISM')
        foldername_ISMIP = [foldername_base '/' henk( i).name];
        break
      end
    end
    % Safety
    if isempty( foldername_ISMIP)
      error(['Couldnt find ISMIP output folder in UFEMISM results folder "' foldername_base '"!'])
    end

    % Find the ISMIP6 output file for total grounded ice mass
    filename_ISMIP = '';
    henk = dir( foldername_ISMIP);
    for i = 1: length( henk)
      if contains( henk( i).name,'limnsw')
        filename_ISMIP = [foldername_ISMIP '/' henk( i).name];
      end
    end
    % Safety
    if isempty( filename_ISMIP)
      error(['Couldnt find limnsw file in ISMIP results folder "' foldername_ISMIP '"!'])
    end

    % Read data
    time2     = ncread( filename_ISMIP,'time') / 360; % 360-day calendar
    VAF_MSLE2 = ncread( filename_ISMIP,'limnsw') / (seawater_density * ocean_area); % expressed in m.s.l.e.

    % DENK DROM - should fix this in the files, first time frame is written twice
    time2     = time2(     2:end);
    VAF_MSLE2 = VAF_MSLE2( 2:end);
    
    if (time2(end)<2300)
      disp([experiment '_' model_version_short ': part 2 runs until ' num2str(time2(end))])
    elseif any(isnan(VAF_MSLE2))
      disp('beep')
    end
    
    %% Combine

    % Save to results structure
    results.(experiment_code).(model_version_short).time         = [time1    ; time2(     2:end)];
    results.(experiment_code).(model_version_short).VAF_MSLE     = [VAF_MSLE1; VAF_MSLE2( 2:end)];
    results.(experiment_code).(model_version_short).VAF_MSLE_rel = results.(experiment_code).(model_version_short).VAF_MSLE - ...
      results.historical.(model_version_short).VAF_MSLE( end);

  end % for mi = 1: length( model_versions)
end % for xi = 1: length( experiments)
  
%% Create plumes
experiment_codes_extra = fields( results);

for xi = 1: length( experiment_codes_extra)
  
  experiment_code = experiment_codes_extra{ xi};
  
  results.(experiment_code).time         = results.(experiment_code).v1.time;
  results.(experiment_code).VAF_MSLE_min = results.(experiment_code).v1.time * 0 + Inf;
  results.(experiment_code).VAF_MSLE_max = results.(experiment_code).v1.time * 0 - Inf;
  results.(experiment_code).VAF_MSLE_av  = results.(experiment_code).v1.time * 0;

  for ti = 1: length( results.(experiment_code).time)

    n = 0;

    for mi = 1: length( model_versions)

      model_version       = model_versions{       mi};
      model_version_short = model_versions_short{ mi};

      tj = find( results.(experiment_code).(model_version_short).time == results.(experiment_code).time( ti));

      if tj > 0
        n = n+1;
        VAF_MSLE = results.(experiment_code).(model_version_short).VAF_MSLE_rel( tj);
        results.(experiment_code).VAF_MSLE_min( ti) = min( results.(experiment_code).VAF_MSLE_min( ti),  VAF_MSLE);
        results.(experiment_code).VAF_MSLE_max( ti) = max( results.(experiment_code).VAF_MSLE_max( ti),  VAF_MSLE);
        results.(experiment_code).VAF_MSLE_av(  ti) =      results.(experiment_code).VAF_MSLE_av(  ti) + VAF_MSLE;
      end

    end % for mi = 1: length( model_versions)
    
    results.(experiment_code).VAF_MSLE_av( ti) = results.(experiment_code).VAF_MSLE_av( ti) / n;
    
  end % for ti = 1: length( results.(experiment_code).time)
end % for xi = 1: length( experiments)

%% Plot

% Empty objects for legend
% colors = [0,0,0; linspecer( length( experiments)-1,'sequential')];
colors = linspecer( length( experiments),'sequential');
for xi = 1: length( experiments)
  line('parent',H.Ax{1,1},'xdata',[],'ydata',[],'linewidth',5,'color',colors( xi,:));
end

% Historical period
color = 'k';

xdata = [results.historical.time        ; flipud( results.historical.time        )];
ydata = [results.historical.VAF_MSLE_min; flipud( results.historical.VAF_MSLE_max)];
patch('parent',H.Ax{1,1},'xdata',xdata,'ydata',-ydata,'facecolor',color,'edgecolor','none','facealpha',0.5);
line( 'parent',H.Ax{1,1},'xdata',results.historical.time,'ydata',-results.historical.VAF_MSLE_min,'color',color,'linewidth',2);
line( 'parent',H.Ax{1,1},'xdata',results.historical.time,'ydata',-results.historical.VAF_MSLE_max,'color',color,'linewidth',2);
line( 'parent',H.Ax{1,1},'xdata',results.historical.time,'ydata',-results.historical.VAF_MSLE_av ,'color',color,'linewidth',3);

% Experiments
for xi = 1: length( experiments)
  
  experiment      = experiments{      xi};
  experiment_code = experiment_codes{ xi};
  
  color = colors( xi,:);
  
  xdata = [results.(experiment_code).time        ; flipud( results.(experiment_code).time        )];
  ydata = [results.(experiment_code).VAF_MSLE_min; flipud( results.(experiment_code).VAF_MSLE_max)];
  patch('parent',H.Ax{1,1},'xdata',xdata,'ydata',-ydata,'facecolor',color,'edgecolor','none','facealpha',0.5);
%   line( 'parent',H.Ax{1,1},'xdata',results.(experiment_code).time,'ydata',-results.(experiment_code).VAF_MSLE_min,'color',color,'linewidth',2);
%   line( 'parent',H.Ax{1,1},'xdata',results.(experiment_code).time,'ydata',-results.(experiment_code).VAF_MSLE_max,'color',color,'linewidth',2);
%   line( 'parent',H.Ax{1,1},'xdata',results.(experiment_code).time,'ydata',-results.(experiment_code).VAF_MSLE_av ,'color',color,'linewidth',3);
  
end
for xi = 1: length( experiments)
  
  experiment      = experiments{      xi};
  experiment_code = experiment_codes{ xi};
  
  color = colors( xi,:);
  
  xdata = [results.(experiment_code).time        ; flipud( results.(experiment_code).time        )];
  ydata = [results.(experiment_code).VAF_MSLE_min; flipud( results.(experiment_code).VAF_MSLE_max)];
%   patch('parent',H.Ax{1,1},'xdata',xdata,'ydata',-ydata,'facecolor',color,'edgecolor','none','facealpha',0.5);
  line( 'parent',H.Ax{1,1},'xdata',results.(experiment_code).time,'ydata',-results.(experiment_code).VAF_MSLE_min,'color',color,'linewidth',2);
  line( 'parent',H.Ax{1,1},'xdata',results.(experiment_code).time,'ydata',-results.(experiment_code).VAF_MSLE_max,'color',color,'linewidth',2);
  line( 'parent',H.Ax{1,1},'xdata',results.(experiment_code).time,'ydata',-results.(experiment_code).VAF_MSLE_av ,'color',color,'linewidth',3);
  
end

% Bars
set( H.Ax{1,2},'xlim',[-0.2,length(experiments)+0.2]);
wb = 0.8;
for xi = 1: length( experiments)
  
  experiment = experiment_codes{ xi};
  experiment_code = experiment_codes{ xi};
  
  color = colors( xi,:);
  
  xlo = (xi-1) + (1-wb)/2;
  xhi = xlo + wb;
  y   = results.(experiment_code).VAF_MSLE_av(  end);
  ylo = results.(experiment_code).VAF_MSLE_min( end);
  yhi = results.(experiment_code).VAF_MSLE_max( end);
  xdata = [xlo,xhi,xhi,xlo];
  ydata = [ylo,ylo,yhi,yhi];
  
  patch('parent',H.Ax{1,2},'xdata',xdata,'ydata',-ydata,'facecolor',color,'edgecolor','none','facealpha',0.5);
  line('parent',H.Ax{1,2},'xdata',[xlo,xhi],'ydata',-[y,y],'color',color,'linewidth',3);
end

% Legend
legend_strings = experiments;
for li = 1: length( legend_strings)
  str = legend_strings{ li};
  i = 1;
  while i < length( str)
    if strcmpi( str( i),'_')
      str = [str( 1:i-1) '\_' str( i+1:end)];
      i = i+2;
    else
      i = i+1;
    end
  end
  legend_strings{ li} = str;
end

legend( H.Ax{1,1},legend_strings,'location','southwest')