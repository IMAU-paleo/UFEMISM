clc
clear all
close all

foldernames = {...
  'MISMIP_mod_DIVA_10km',...
  'MISMIP_mod_hybrid_10km',...
  'MISMIP_mod_DIVA_sans_10km',...
  'MISMIP_mod_hybrid_sans_10km'};

%% Read model output

% for fi = 1: length( foldernames)
%   
%   foldername = foldernames{ fi};
%   
%   timeframes = get_UFEMISM_filelist( foldername, 'ANT');
%   ntf = length( timeframes);
%   
%   results( fi).time   = zeros( ntf,1);
%   results( fi).xGL_av = zeros( ntf,1);
% 
%   filename_prev = '';
%   
%   for tfi = 1: ntf
%   
%     % Get info for this timeframe
%     time                 = timeframes( tfi).time;
%     ti                   = timeframes( tfi).ti;
%     filename_restart     = timeframes( tfi).filename_restart;
%     filename_help_fields = timeframes( tfi).filename_help_fields;
% 
%     % If necessary, read the new mesh
%     if ~strcmpi( filename_prev, filename_restart)
%       filename_prev = filename_restart;
%       mesh = read_mesh_from_file( filename_restart);
%     end
% 
%     % Read data fields for this timeframe
%     Hi = ncread( filename_restart,'Hi',[1,ti],[Inf,1]);
%     Hb = ncread( filename_restart,'Hb',[1,ti],[Inf,1]);
%     SL = ncread( filename_restart,'SL',[1,ti],[Inf,1]);
%     
%     % Calculate thickness above flotation
%     ice_density      = 910;
%     seawater_density = 1028;
%     TAF = Hi - max(0, (SL - Hb) * (seawater_density / ice_density));
%     
%     % Calculate GL as TAF=0 contour
%     C_GL = mesh_contour( mesh,TAF,0);
%     
%     % Calculate average GL position
%     xGL = sqrt( C_GL(:,1).^2 + C_GL(:,2).^2);
%     xGL_av = mean( xGL);
%     
%     % Save in results struct
%     results( fi).time(   tfi) = time;
%     results( fi).xGL_av( tfi) = xGL_av;
%     
%   end
% end
% 
% save('tempdata_withvssans.mat','results');
load('tempdata_withvssans.mat');

%% Plot results

wa = 600;
ha = 400;

margin_left   = 110;
margin_right  = 25;
margin_bottom = 80;
margin_top    = 25;

wf = margin_left + wa + margin_right;
hf = margin_bottom + ha + margin_top;

H.Fig = figure('color','w','position',[300,300,wf,hf]);
H.Ax1 = axes('parent',H.Fig,'units','pixels','position',[margin_left,margin_bottom,wa,ha],...
  'fontsize',24,'xgrid','on','ygrid','on','xlim',[0,45],'ylim',[750,1050]);

xlabel(H.Ax1,'Time (kyr)')
ylabel(H.Ax1,'x_{GL} (km)')

% Empty line objects for legend
c_with = [1.0,0.2,0.2];
c_sans = [1.0,0.7,0.4];

line(H.Ax1,'xdata',[],'ydata',[],'color',c_with,'linewidth',3,'linestyle','-');   % DIVA
line(H.Ax1,'xdata',[],'ydata',[],'color',c_with,'linewidth',3,'linestyle','--');  % Hybrid
line(H.Ax1,'xdata',[],'ydata',[],'color',c_sans,'linewidth',3,'linestyle','-');   % DIVA sans
line(H.Ax1,'xdata',[],'ydata',[],'color',c_sans,'linewidth',3,'linestyle','--');  % Hybrid sans

% Plot results
line('parent',H.Ax1,'xdata',results( 1).time/1e3,'ydata',results( 1).xGL_av / 1e3,'linewidth',3,'color',c_with,'linestyle','-');  % DIVA
line('parent',H.Ax1,'xdata',results( 2).time/1e3,'ydata',results( 2).xGL_av / 1e3,'linewidth',3,'color',c_with,'linestyle','--'); % Hybrid
line('parent',H.Ax1,'xdata',results( 3).time/1e3,'ydata',results( 3).xGL_av / 1e3,'linewidth',3,'color',c_sans,'linestyle','-');  % DIVA sans
line('parent',H.Ax1,'xdata',results( 4).time/1e3,'ydata',results( 4).xGL_av / 1e3,'linewidth',3,'color',c_sans,'linestyle','--'); % Hybrid sans

legend(H.Ax1,'DIVA','hybrid','DIVA sans','hybrid sans','location','northwest')