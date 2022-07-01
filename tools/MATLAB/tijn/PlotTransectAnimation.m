clc
clear all
close all

foldername = '../results_20201111_001';

ylim = [-2000,6000];

tstart = -25000;
tstop  = 0;
dt     = 100;

%% Create graphics objects

ah = 400;
aw = 700;

fh = 100+ah+75;
fw = 100+aw+25;

H.Fig = figure('position',[200 200 fw fh],'color','w');
H.Ax  = axes('units','pixels','position',[100 100 aw ah],'box','on','fontsize',28,...
  'ylim',ylim,'color',[0.9,0.95,1.0],'ytick',-10000:1000:10000,'xgrid','on','ygrid','on');

H.patch_bed = patch('vertices',[],'faces',[],'facecolor',[0.6,0.3,0.1],'edgecolor','none');
H.patch_ice = patch('vertices',[],'faces',[],'facecolor','w','edgecolor','none');
H.line_ice_top = line('xdata',[],'ydata',[]);
H.line_ice_bot = line('xdata',[],'ydata',[]);

%% Plot all output files for this region
nfile = 1;
filename = [foldername '/results_ANT_00001.nc'];
while exist(filename,'file')
  
  % Update transect x data
  trans.x = ReadTransectX( filename);
  
  set(H.Ax,'xlim',[min(trans.x),max(trans.x)]);
  
  % Go through all this file's time frames
  time = ncread(filename,'time');
  
  if ~isempty(time) && max(time)>=tstart && min(time)<=tstop
    for ti = 1:1:length(time)

      if mod(time(ti),dt)>0; continue; end
      if time(ti)<tstart || time(ti)>tstop; continue; end

      title(H.Ax,['Time: ' num2str(round(time(ti))) ' y'])

      % Read data
      trans.Hi = ReadTransectData( filename, 'Hi', ti);
      trans.Hb = ReadTransectData( filename, 'Hb', ti);
      trans.Hs = ReadTransectData( filename, 'Hs', ti);

      % Update patches
      % Bed
      xdata = [trans.x;flipud(trans.x)];
      ydata = [zeros(size(trans.x))+ylim(1);flipud(trans.Hb)];
      set(H.patch_bed,'xdata',xdata,'ydata',ydata);

      % Ice
      xdata = [trans.x;flipud(trans.x)];
      ice_bot = trans.Hs - trans.Hi;
      ice_top = trans.Hs;
      ydata = [ice_bot;flipud(ice_top)];
      set(H.patch_ice,'xdata',xdata,'ydata',ydata);
      set(H.line_ice_top,'xdata',trans.x,'ydata',ice_top);
      set(H.line_ice_bot,'xdata',trans.x,'ydata',ice_bot);

      drawnow('update');
%       pause

    end
  end
  
  nfile = nfile+1;
  if     nfile<10
    nfilestr = ['0000' num2str(nfile)];
  elseif nfile<100
    nfilestr = ['000' num2str(nfile)];
  elseif nfile<1000
    nfilestr = ['00' num2str(nfile)];
  elseif nfile<10000
    nfilestr = ['0' num2str(nfile)];
  else
    nfilestr = num2str(nfile);
  end
filename = [foldername '/results_ANT_' nfilestr '.nc'];
end