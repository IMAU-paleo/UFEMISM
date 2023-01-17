clc
clear all
close all

%% Read other model data
foldername ='/Users/berends/Documents/Datasets/ISMIP-HOM/tc-2-95-2008-supplement/ismip_all/';

model_types = {...
  'aas1','FS';...
  'aas2','FS';...
  'ahu1','HO';...
  'ahu2','HO';...
  'bds1','HO';...
  'cma1','FS';...
  'cma2','HO';...
  'dpo1','HO';...
  'fpa1','HO';...
  'fpa2','FS';...
  'fsa1','HO';...
  'ghg1','FS';...
  'jvj1','FS';...
  'lpe1','HO';...
  'mbr1','HO';...
  'mmr1','FS';...
  'mtk1','HO';...
  'oga1','FS';...
  'oso1','HO';...
  'rhi1','FS';...
  'rhi2','HO';...
  'rhi3','FS';...
  'rhi4','HO';...
  'rhi5','HO';...
  'spr1','FS';...
  'ssu1','FS';...
  'tpa1','HO';...
  'yko1','FS'};

models = dir(foldername);
models = models(3:end);

% experiments.a160 = [];
% experiments.a080 = [];
% experiments.a040 = [];
% experiments.a020 = [];
% experiments.a010 = [];
% experiments.a005 = [];
% 
% for mi = 1: length(models)
%   
%   modeldata = dir([foldername models(mi).name]);
%   modeldata = modeldata(3:end);
%   
%   % Go over all experiments, check if this model has them.
%   flds = fields(experiments);
%   for xi = 1:length(flds)
%     ex = flds{xi};
%     
%     for di = 1:length(modeldata)
%       mdname = modeldata(di).name;
%       str = [ex '.txt'];
%       if length(mdname) >= length(str)
%         if strcmpi(mdname(end-length(str)+1:end),str)
%           % This is the experiment from this model
%           
%           disp(['Reading data from model ' models(mi).name ', experiment ' ex])
%           
%           fid = fopen([foldername models(mi).name '/' mdname]);
%           temp = textscan(fid,'%s','delimiter','\n','MultipleDelimsAsOne',1); temp = temp{1};
%           fclose(fid);
%           
%           n = length(temp);
%           nx = sqrt(n);
%           if nx-floor(nx)>0
%             error('whaa!')
%           end
%           x_vec = zeros(n,1);
%           y_vec = zeros(n,1);
%           u_vec = zeros(n,1);
%           v_vec = zeros(n,1);
%           w_vec = zeros(n,1);
%           for i = 1:n
%             temp2 = textscan(temp{i},'%f %f %f %f %f %f %f %f');
%             x_vec(i) = temp2{1};
%             y_vec(i) = temp2{2};
%             u_vec(i) = temp2{3};
%             v_vec(i) = temp2{4};
%             w_vec(i) = temp2{5};
%           end
%           
%           experiments.(ex).(models(mi).name).x = reshape(x_vec,[nx,nx]);
%           experiments.(ex).(models(mi).name).y = reshape(y_vec,[nx,nx]);
%           experiments.(ex).(models(mi).name).u = reshape(u_vec,[nx,nx]);
%           experiments.(ex).(models(mi).name).v = reshape(v_vec,[nx,nx]);
%           experiments.(ex).(models(mi).name).w = reshape(w_vec,[nx,nx]);
%           
%           if (experiments.(ex).(models(mi).name).x(1,1) == experiments.(ex).(models(mi).name).x(end,1))
%             experiments.(ex).(models(mi).name).x = experiments.(ex).(models(mi).name).x';
%             experiments.(ex).(models(mi).name).y = experiments.(ex).(models(mi).name).y';
%             experiments.(ex).(models(mi).name).u = experiments.(ex).(models(mi).name).u';
%             experiments.(ex).(models(mi).name).v = experiments.(ex).(models(mi).name).v';
%             experiments.(ex).(models(mi).name).w = experiments.(ex).(models(mi).name).w';
%           end
%           
%         end
%       end
%     end
%   end
%   
% end
% 
% save('ISMIP_HOM_A_ensemble.mat','experiments');
load('ISMIP_HOM_A_ensemble.mat')

%% Plot Pattyn2008 ensemble

wa = 600;
ha = 400;
margin_left  = 100;
margin_right = 25;
margin_bot   = 80;
margin_top   = 25;
wf = margin_left + wa + margin_right;
hf = margin_bot  + ha + margin_top;

H.Fig = figure('position',[200,200,wf,hf],'color','w');

H.Ax = axes('units','pixels','position',[margin_left,margin_bot,wa,ha],...
  'fontsize',24,'xlim',[0,1],'xgrid','on','ygrid','on');

xlabel(H.Ax,'x / L')
ylabel(H.Ax,'Velocity (m/yr)')

set(H.Ax,'ylim',[0,125]);

% Legend
c_UFE = [1.0,0.2,0.2];
c_FS     = [0.2,0.5,1.0];
c_HO     = [0.1,0.7,0.3];
line(H.Ax(1),'xdata',[],'ydata',[],'color',c_UFE,'linewidth',3,'linestyle','-');
patch(H.Ax(1),'vertices',[],'faces',[],'facecolor',c_FS,'edgecolor','none','facealpha',0.7)
patch(H.Ax(1),'vertices',[],'faces',[],'facecolor',c_HO,'edgecolor','none','facealpha',0.7)
line(H.Ax(1),'xdata',[],'ydata',[],'color',c_FS,'linewidth',3);
line(H.Ax(1),'xdata',[],'ydata',[],'color',c_HO,'linewidth',3);

% Results
for a = 1
  
  if a == 1
    ex = 'a160';
  elseif a == 2
    ex = 'a080';
  elseif a == 3
    ex = 'a040';
  elseif a == 4
    ex = 'a020';
  elseif a == 5
    ex = 'a010';
  elseif a == 6
    ex = 'a005';
  end
  
  x_FS = linspace(0,1,200)';
  u_FS = [];
  u_HO = [];
  
  patch_HO = patch(H.Ax(a),'xdata',[],'ydata',[],'facecolor',c_HO,'edgecolor','none','facealpha',0.5);
  patch_FS = patch(H.Ax(a),'xdata',[],'ydata',[],'facecolor',c_FS,'edgecolor','none','facealpha',0.7);
  line_HO  = line('parent',H.Ax(a),'xdata',[],'ydata',[],'color',c_HO,'linewidth',3);
  line_FS  = line('parent',H.Ax(a),'xdata',[],'ydata',[],'color',c_FS,'linewidth',3);
  
  flds = fields(experiments.(ex));
  for mi = 1:length(flds)
    m = flds{mi};
    
    x = experiments.(ex).(m).x(:,1);
    y = experiments.(ex).(m).y(1,:)';
    yi = round(0.25*length(y));
    u = experiments.(ex).(m).u(:,yi);

    % Determine if this model is FS or HO
    FS = false;
    HO = false;
    for mii = 1:size(model_types,1)
      if strcmpi(model_types{mii,1},m)
        if strcmpi(model_types{mii,2},'FS')
          FS = true;
        else
          HO = true;
        end
      end
    end
    if ~(FS || HO)
      for mii = 1:size(model_types,1)
        if strcmpi(model_types{mii,1}(1:3),m(1:3))
          if strcmpi(model_types{mii,2},'FS')
            FS = true;
          else
            HO = true;
          end
        end
      end
    end
    if ~(FS || HO)
      % Unknown model?
      continue
    end

    % Add to data ranges for HO/FS models
    up = interp1(x,u,x_FS);
    if FS
      u_FS(:,end+1) = up;
    else
      u_HO(:,end+1) = up;
    end
    
  end
  
  % Missing data points
  m = true(size(x_FS));
  for i = 1:length(x_FS)
    if sum(isnan(u_FS(i,:)))+sum(isnan(u_HO(i,:)))>0
      m(i) = false;
    end
  end
  u_FS = u_FS(m,:);
  u_HO = u_HO(m,:);
  x_FS = x_FS(m);
  
  % ISMIP-HOM ensemble data
  uav_FS = mean(u_FS,2);
  uav_HO = mean(u_HO,2);
  sigma_FS = zeros(size(x_FS));
  sigma_HO = zeros(size(x_FS));
  for i = 1:size(u_FS,1)
    sigma_FS(i) = std(u_FS(i,:));
    sigma_HO(i) = std(u_HO(i,:));
    if isnan(sigma_FS(i)); sigma_FS(i) = 0; end
    if isnan(sigma_HO(i)); sigma_HO(i) = 0; end
  end
  umin_FS = uav_FS - sigma_FS;
  umax_FS = uav_FS + sigma_FS;
  umin_HO = uav_HO - sigma_HO;
  umax_HO = uav_HO + sigma_HO;
  
  xdata = [x_FS;flipud(x_FS)];
  ydata = [umin_FS;flipud(umax_FS)];
  set(patch_FS,'xdata',xdata,'ydata',ydata)
  
  xdata = [x_FS;flipud(x_FS)];
  ydata = [umin_HO;flipud(umax_HO)];
  set(patch_HO,'xdata',xdata,'ydata',ydata)
  
  set(line_HO,'xdata',x_FS,'ydata',uav_HO)
  set(line_FS,'xdata',x_FS,'ydata',uav_FS)
  
end

%% Read and plot Matlab-UFEMISM DIVA/BPA results

nzs   = [8,12,16];
Rstrs = {'0p1','0p2','0p5','01','02','05','10'};

for nzi = 1: length( nzs)
  nz = nzs( nzi);
  for ri = 1: length(Rstrs)

    L = 160e3;

    trans.n  = 100;
    trans.xt = linspace( -L/2, L/2, trans.n)';
    trans.yt = zeros( trans.n,1) + L/4;
    
    foldername = ['ISMIP_HOM_A_BPA_nz' num2str(nz) '_R' Rstrs{ri} '_160'];
    filename = [foldername '/help_fields_ANT_00001.nc'];
    
    if ~exist(filename,'file'); continue; end
    time = ncread( filename,'time');
    if isempty( time); continue; end
    
    mesh     = read_mesh_from_file( filename);
    trans.A  = calc_transect_matrix_a( mesh, trans.xt, trans.yt);
    time     = ncread( filename,'time');
    ti       = length( time);
    u_surf_a = ncread( filename,'u_3D',[1,1,ti],[Inf,1,1]);
    trans.u_surf = trans.A * u_surf_a;
    line('parent',H.Ax,'xdata',fliplr(linspace(0,1,trans.n)),'ydata',trans.u_surf,'linewidth',1,...
      'color',c_UFE,'linestyle','-');

  end
end

legend(H.Ax,'BPA','Full-Stokes','Higher-Order','FS mean','HO mean','location','northwest')