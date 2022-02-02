clc
clear all
close all

foldernames = {...
  'MISMIP_mod_hybrid_64km',...
  'MISMIP_mod_hybrid_40km',...
  'MISMIP_mod_hybrid_32km',...
  'MISMIP_mod_hybrid_20km',...
  'MISMIP_mod_hybrid_16km',...
  'MISMIP_mod_hybrid_10km'};

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
% save('tempdat259768a_hybrid.mat','results');
load('tempdata_hybrid.mat');

%% Plot results

wa1 = 600;
wa2 = 400;
ha = 400;

margin_left   = 110;
margin_mid    = 100;
margin_right  = 25;
margin_bottom = 80;
margin_top    = 25;

wf = margin_left + wa1 + margin_mid + wa2 + margin_right;
hf = margin_bottom + ha + margin_top;

H.Fig = figure('color','w','position',[300,300,wf,hf]);
H.Ax1 = axes('parent',H.Fig,'units','pixels','position',[margin_left,margin_bottom,wa1,ha],...
  'fontsize',24,'xgrid','on','ygrid','on','xlim',[0,45],'ylim',[750,1050]);

%% Timeseries
colors = flipud( parula( length(foldernames)));

% Empty line objects for legend
for fi = 1: length( foldernames)
  line('parent',H.Ax1,'xdata',[],'ydata',[],'color',colors( fi,:),'linewidth',3);
end

% Plot results
for fi = 1: length( foldernames)
  line('parent',H.Ax1,'xdata',results( fi).time/1e3,'ydata',results( fi).xGL_av / 1e3,...
    'color',colors( fi,:),'linewidth',3);
end

xlabel(H.Ax1,'Time (kyr)')
ylabel(H.Ax1,'x_{GL} (km)')

% legend(H.Ax1,'64 km','40 km','32 km','20 km','16 km','10 km','8 km','5 km','4 km','location','northwest')
legend(H.Ax1,'64 km','40 km','32 km','20 km','16 km','10 km','location','northwest')

%% Hysteresis

H.Ax2 = axes('parent',H.Fig,'units','pixels','position',[margin_left+wa1+margin_mid,margin_bottom,wa2,ha],...
  'fontsize',24,'xgrid','on','ygrid','on','xscale','log','yscale','log',...
  'xlim',[2,128],'xtick',[4,8,16,32,64],'ylim',[2,128],'ytick',[4,8,16,32,64]);

xlabel(H.Ax2,'Resolution (km)')
ylabel(H.Ax2,'GL hysteresis (km)')

for fi = 1: length(foldernames)
  
  ti1 = find( results( fi).time == 15000);
  ti2 = find( results( fi).time == 45000);
  
  xGL1 = results( fi).xGL_av( ti1);
  xGL2 = results( fi).xGL_av( ti2);
  
  dxGL( fi) = abs( xGL2 - xGL1) / 1e3;
  
end

resolutions = [64,40,32,20,16,10];%,8,5,4];

line('parent',H.Ax2,'xdata',[],'ydata',[],'color','b','linewidth',2)

line('parent',H.Ax2,'xdata',[2,128],'ydata',[2,128],'linestyle','--')
line('parent',H.Ax2,'xdata',resolutions,'ydata',dxGL,'linestyle','none','marker','o','markerfacecolor','b','markeredgecolor','b','markersize',8);

% Loglinear fit
p = polyfit( log(resolutions), log(dxGL), 1);
dxGL_fit = exp(polyval(p,log([2,128])));
line('parent',H.Ax2,'xdata',[2,128],'ydata',dxGL_fit,'color','b','linestyle','-','linewidth',2)

legend(['O( R^{' num2str(round(p(1)*100)/100) '})'],'location','northwest')