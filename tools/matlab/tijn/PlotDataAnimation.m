clc
clear all
close all

foldername = '../results_20210316_001';
regionname = 'ANT';

dataname = 'Hi';

cmap = parula(256);
% clim = [0 1e6];
edgecolor = 'w';

tstart = -60000;
tstop  = 60000;
dt     = 5000;  % Time between output frames (years)

%% Create graphics objects

V   = ncread([foldername '/results_' regionname '_00001.nc'],'V');
Tri = ncread([foldername '/results_' regionname '_00001.nc'],'Tri');

xr = max(V(:,1)) - min(V(:,1));
yr = max(V(:,2)) - min(V(:,2));

ah = 500;
aw = ah*xr/yr;

fh = ah+100;
fw = aw+150;

H.Fig = figure('position',[900 450 fw fh],'color','w');
H.Ax  = axes('units','pixels','position',[25 25 aw ah],'xtick',[],'ytick',[],'box','on',...
  'xlim',[min(V(:,1)) max(V(:,1))],'ylim',[min(V(:,2)) max(V(:,2))],'fontsize',28);
if exist('clim','var'); set(H.Ax,'clim',clim); end

H.Patch = patch('parent',H.Ax,'vertices',V,'faces',Tri,'facecolor','none','edgecolor',edgecolor);

colormap(H.Ax,cmap);

H.Cbar = colorbar(H.Ax,'location','eastoutside');
set(H.Ax,'position',[25 25 aw ah]);
set(H.Ax,'units','normalized');

%% Plot all output files for this region
nfile = 1;
filename = [foldername '/results_' regionname '_00001.nc'];
while exist(filename,'file')
  
  % Update mesh data
  V   = ncread(filename,'V');
  Tri = ncread(filename,'Tri');  
  
  set(H.Patch,'vertices',V,'faces',Tri);
  
  % Go through all this file's time frames
  time = ncread(filename,'time');
  
  if ~isempty(time) && max(time)>=tstart && min(time)<=tstop
    for ti = 1:1:length(time)

      if mod(time(ti),dt)>0; continue; end
      if time(ti)<tstart || time(ti)>tstop; continue; end

      title(H.Ax,['Time: ' num2str(round(time(ti))) ' y'])

      % Update model data
      f = ncinfo(filename,dataname);
      if length(f.Size)==1
        d = ncread(filename, dataname, 1, Inf); 
      elseif length(f.Size)==2
        d = ncread(filename, dataname, [1 ti],[Inf 1]);    
      else
        d = ncread(filename, dataname, [1 f.Size(2) ti],[Inf 1 1]); d = mean(d,2);
      end

      if strcmpi(dataname,'R')
        d = -log(d);
      end

%       d1 = ncread(filename,'U_SSA',[1,ti],[Inf,1]);
%       d2 = ncread(filename,'V_SSA',[1,ti],[Inf,1]);
%       d = log(sqrt(d1.^2+d2.^2));

      set(H.Patch,'facevertexcdata',d,'facecolor','interp');
      drawnow('update');
      pause

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
filename = [foldername '/results_' regionname '_' nfilestr '.nc'];
end