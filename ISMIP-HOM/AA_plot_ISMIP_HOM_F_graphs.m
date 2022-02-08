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

% experiments.f000 = [];
% experiments.f001 = [];
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
%           x_vec  = zeros(n,1);
%           y_vec  = zeros(n,1);
%           Hs_vec = zeros(n,1);
%           u_vec  = zeros(n,1);
%           v_vec  = zeros(n,1);
%           w_vec  = zeros(n,1);
%           for i = 1:n
%             temp2 = textscan(temp{i},'%f %f %f %f %f %f');
%             x_vec( i) = temp2{1};
%             y_vec( i) = temp2{2};
%             Hs_vec(i) = temp2{3};
%             u_vec( i) = temp2{4};
%             v_vec( i) = temp2{5};
%             w_vec( i) = temp2{6};
%           end
%           
%           x_vec = (x_vec - min(x_vec)) / range(x_vec);
%           
%           experiments.(ex).(models(mi).name).x  = reshape(x_vec, [nx,nx]);
%           experiments.(ex).(models(mi).name).y  = reshape(y_vec, [nx,nx]);
%           experiments.(ex).(models(mi).name).Hs = reshape(Hs_vec,[nx,nx]);
%           experiments.(ex).(models(mi).name).u  = reshape(u_vec, [nx,nx]);
%           experiments.(ex).(models(mi).name).v  = reshape(v_vec, [nx,nx]);
%           experiments.(ex).(models(mi).name).w  = reshape(w_vec, [nx,nx]);
%           
%           if (experiments.(ex).(models(mi).name).x(1,1) == experiments.(ex).(models(mi).name).x(end,1))
%             experiments.(ex).(models(mi).name).x  = experiments.(ex).(models(mi).name).x';
%             experiments.(ex).(models(mi).name).y  = experiments.(ex).(models(mi).name).y';
%             experiments.(ex).(models(mi).name).Hs = experiments.(ex).(models(mi).name).Hs';
%             experiments.(ex).(models(mi).name).u  = experiments.(ex).(models(mi).name).u';
%             experiments.(ex).(models(mi).name).v  = experiments.(ex).(models(mi).name).v';
%             experiments.(ex).(models(mi).name).w  = experiments.(ex).(models(mi).name).w';
%           end
%           
%         end
%       end
%     end
%   end
%   
% end
% 
% save('ISMIP_HOM_F_ensemble.mat','experiments')
load('ISMIP_HOM_F_ensemble.mat')

%% Read my own results

% DIVA
UFEMISM.F0.DIVA.foldername        = '/Users/berends/Documents/Models/UFEMISM/ISMIP-HOM/ISMIP_HOM_F_0_DIVA';
UFEMISM.F0.hybrid.foldername      = '/Users/berends/Documents/Models/UFEMISM/ISMIP-HOM/ISMIP_HOM_F_0_hybrid';
UFEMISM.F0.DIVA_sans.foldername   = '/Users/berends/Documents/Models/UFEMISM/ISMIP-HOM/ISMIP_HOM_F_0_DIVA_sans';
UFEMISM.F0.hybrid_sans.foldername = '/Users/berends/Documents/Models/UFEMISM/ISMIP-HOM/ISMIP_HOM_F_0_hybrid_sans';

UFEMISM.F1.DIVA.foldername        = '/Users/berends/Documents/Models/UFEMISM/ISMIP-HOM/ISMIP_HOM_F_1_DIVA';
UFEMISM.F1.hybrid.foldername      = '/Users/berends/Documents/Models/UFEMISM/ISMIP-HOM/ISMIP_HOM_F_1_hybrid';
UFEMISM.F1.DIVA_sans.foldername   = '/Users/berends/Documents/Models/UFEMISM/ISMIP-HOM/ISMIP_HOM_F_1_DIVA_sans';
UFEMISM.F1.hybrid_sans.foldername = '/Users/berends/Documents/Models/UFEMISM/ISMIP-HOM/ISMIP_HOM_F_1_hybrid_sans';

exps = {'F0','F1'};
icedyns = {'DIVA','hybrid','DIVA_sans','hybrid_sans'};
for xi = 1: length( exps)
  for icei = 1: length( icedyns)
    
    foldername           = UFEMISM.(exps{xi}).(icedyns{icei}).foldername;
    filename_restart     = [foldername '/restart_grid_ANT.nc'];
    filename_help_fields = [foldername '/help_fields_grid_ANT.nc'];
    
    time = ncread( filename_restart,'time');
    ti = length( time);

    x  = ncread( filename_restart,'x');
    y  = ncread( filename_restart,'y');
    j = find(y==0);
    
    Hs = ncread( filename_restart    ,'Hs'  ,[1,j,  ti],[Inf,1,  1]) - ncread( filename_restart    ,'Hs'  ,[1,j,  1],[Inf,1,  1]);
    u  = ncread( filename_help_fields,'u_3D',[1,j,1,ti],[Inf,1,1,1]);
    
    % Scale and crop x
    x = 3*(x - min(x)) / range(x) - 1;
    Hs(x>1) = [];
    u( x>1) = [];
    x( x>1) = [];
    Hs(x<0) = [];
    u( x<0) = [];
    x( x<0) = [];

    UFEMISM.(exps{xi}).(icedyns{icei}).x  = x;
    UFEMISM.(exps{xi}).(icedyns{icei}).Hs = Hs;
    UFEMISM.(exps{xi}).(icedyns{icei}).u  = u;
    
  end
end

%% Plot results
wa = 400;
ha = 300;
margin_left   = 90;
margin_right  = 90;
margin_mid_x  = 25;
margin_bot    = 100;
margin_mid_y  = 50;
margin_top    = 50;
wf = margin_left + wa + margin_mid_x + wa + margin_right;
hf = margin_bot + ha + margin_mid_y + ha + margin_top;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
x1 = margin_left;
x2 = margin_left + wa + margin_mid_x;
y1 = margin_bot + ha + margin_mid_y;
y2 = margin_bot;

H.Ax(1) = axes('units','pixels','position',[x1,y1,wa,ha],'xticklabels',[]);
H.Ax(2) = axes('units','pixels','position',[x2,y1,wa,ha],'xticklabels',[],'yaxislocation','right');
H.Ax(3) = axes('units','pixels','position',[x1,y2,wa,ha]);
H.Ax(4) = axes('units','pixels','position',[x2,y2,wa,ha],'yaxislocation','right');

for a = 1:4
  set(H.Ax(a),'fontsize',24,'xlim',[0,1],'xgrid','on','ygrid','on');
end

title(H.Ax(1),'No-slip bed')
title(H.Ax(2),'No-slip bed')
title(H.Ax(3),'Slip bed')
title(H.Ax(4),'Slip bed')

ylabel(H.Ax(1),'Surface elevation (m)')
ylabel(H.Ax(3),'Surface elevation (m)')
ylabel(H.Ax(2),'Velocity (m/yr)')
ylabel(H.Ax(4),'Velocity (m/yr)')

xlabel(H.Ax(3),'x / L')
xlabel(H.Ax(4),'x / L')

% Legend
c_with = [1.0,0.2,0.2];
c_sans = [1.0,0.7,0.4];
c_FS   = [0.2,0.5,1.0];
c_HO   = [0.1,0.7,0.3];
line(H.Ax(1),'xdata',[],'ydata',[],'color',c_with,'linewidth',3,'linestyle','-');   % DIVA
line(H.Ax(1),'xdata',[],'ydata',[],'color',c_with,'linewidth',3,'linestyle','--');  % Hybrid
line(H.Ax(1),'xdata',[],'ydata',[],'color',c_sans,'linewidth',3,'linestyle','-');   % DIVA sans
line(H.Ax(1),'xdata',[],'ydata',[],'color',c_sans,'linewidth',3,'linestyle','--');  % Hybrid sans
patch(H.Ax(1),'vertices',[],'faces',[],'facecolor',c_FS,'edgecolor','none','facealpha',0.7)
patch(H.Ax(1),'vertices',[],'faces',[],'facecolor',c_HO,'edgecolor','none','facealpha',0.7)
line(H.Ax(1),'xdata',[],'ydata',[],'color',c_FS,'linewidth',3);
line(H.Ax(1),'xdata',[],'ydata',[],'color',c_HO,'linewidth',3);

% Results
for a = 1:4
  
  if a == 1
    ex = 'f000';
  elseif a == 2
    ex = 'f000';
  elseif a == 3
    ex = 'f001';
  elseif a == 4
    ex = 'f001';
  end
  
  x_FS = linspace(0,1,200)';
  u_FS = [];
  u_HO = [];
  Hs_FS = [];
  Hs_HO = [];
  
  patch_HO = patch(H.Ax(a),'xdata',[],'ydata',[],'facecolor',c_HO,'edgecolor','none','facealpha',0.5);
  patch_FS = patch(H.Ax(a),'xdata',[],'ydata',[],'facecolor',c_FS,'edgecolor','none','facealpha',0.7);
  line_HO  = line('parent',H.Ax(a),'xdata',[],'ydata',[],'color',c_HO,'linewidth',3);
  line_FS  = line('parent',H.Ax(a),'xdata',[],'ydata',[],'color',c_FS,'linewidth',3);
  
  flds = fields(experiments.(ex));
  for mi = 1:length(flds)
    m = flds{mi};
    
    x = experiments.(ex).(m).x(:,1);
    y = experiments.(ex).(m).y(1,:)';
    yi = round(0.5*length(y));
    u = experiments.(ex).(m).u(:,yi);
    Hs = experiments.(ex).(m).Hs(:,yi);

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
    Hsp = interp1(x,Hs,x_FS);
    if FS
      u_FS(:,end+1) = up;
      Hs_FS(:,end+1) = Hsp;
    else
      u_HO(:,end+1) = up;
      Hs_HO(:,end+1) = Hsp;
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
  Hs_FS = Hs_FS(m,:);
  Hs_HO = Hs_HO(m,:);
  x_FS = x_FS(m);
  
  % ISMIP-HOM ensemble data
  uav_FS = mean(u_FS,2);
  uav_HO = mean(u_HO,2);
  Hsav_FS = mean(Hs_FS,2);
  Hsav_HO = mean(Hs_HO,2);
  sigmau_FS = zeros(size(x_FS));
  sigmau_HO = zeros(size(x_FS));
  sigmaHs_FS = zeros(size(x_FS));
  sigmaHs_HO = zeros(size(x_FS));
  for i = 1:size(u_FS,1)
    sigmau_FS(i) = std(u_FS(i,:));
    sigmau_HO(i) = std(u_HO(i,:));
    sigmaHs_FS(i) = std(Hs_FS(i,:));
    sigmaHs_HO(i) = std(Hs_HO(i,:));
    if isnan(sigmau_FS(i)); sigmau_FS(i) = 0; end
    if isnan(sigmau_HO(i)); sigmau_HO(i) = 0; end
    if isnan(sigmaHs_FS(i)); sigmaHs_FS(i) = 0; end
    if isnan(sigmaHs_HO(i)); sigmaHs_HO(i) = 0; end
  end
  umin_FS = uav_FS - sigmau_FS;
  umax_FS = uav_FS + sigmau_FS;
  umin_HO = uav_HO - sigmau_HO;
  umax_HO = uav_HO + sigmau_HO;
  Hsmin_FS = Hsav_FS - sigmaHs_FS;
  Hsmax_FS = Hsav_FS + sigmaHs_FS;
  Hsmin_HO = Hsav_HO - sigmaHs_HO;
  Hsmax_HO = Hsav_HO + sigmaHs_HO;
  
  if (a == 2 || a == 4)
  
    xdata = [x_FS;flipud(x_FS)];
    ydata = [umin_FS;flipud(umax_FS)];
    set(patch_FS,'xdata',xdata,'ydata',ydata)

    xdata = [x_FS;flipud(x_FS)];
    ydata = [umin_HO;flipud(umax_HO)];
    set(patch_HO,'xdata',xdata,'ydata',ydata)

    set(line_HO,'xdata',x_FS,'ydata',uav_HO)
    set(line_FS,'xdata',x_FS,'ydata',uav_FS)
  
  else
  
    xdata = [x_FS;flipud(x_FS)];
    ydata = [Hsmin_FS;flipud(Hsmax_FS)];
    set(patch_FS,'xdata',xdata,'ydata',ydata)

    xdata = [x_FS;flipud(x_FS)];
    ydata = [Hsmin_HO;flipud(Hsmax_HO)];
    set(patch_HO,'xdata',xdata,'ydata',ydata)

    set(line_HO,'xdata',x_FS,'ydata',Hsav_HO)
    set(line_FS,'xdata',x_FS,'ydata',Hsav_FS)
  
  end
  
end

%% Plot UFEMISM results
line('parent',H.Ax(1),'xdata',UFEMISM.F0.DIVA.x       ,'ydata',UFEMISM.F0.DIVA.Hs       ,'linewidth',3,'color',c_with,'linestyle','-');  % DIVA
line('parent',H.Ax(1),'xdata',UFEMISM.F0.hybrid.x     ,'ydata',UFEMISM.F0.hybrid.Hs     ,'linewidth',3,'color',c_with,'linestyle','--'); % hybrid
line('parent',H.Ax(1),'xdata',UFEMISM.F0.DIVA_sans.x  ,'ydata',UFEMISM.F0.DIVA_sans.Hs  ,'linewidth',3,'color',c_sans,'linestyle','-');  % DIVA_sans
line('parent',H.Ax(1),'xdata',UFEMISM.F0.hybrid_sans.x,'ydata',UFEMISM.F0.hybrid_sans.Hs,'linewidth',3,'color',c_sans,'linestyle','--'); % hybrid_sans

line('parent',H.Ax(2),'xdata',UFEMISM.F0.DIVA.x       ,'ydata',UFEMISM.F0.DIVA.u        ,'linewidth',3,'color',c_with,'linestyle','-');  % DIVA
line('parent',H.Ax(2),'xdata',UFEMISM.F0.hybrid.x     ,'ydata',UFEMISM.F0.hybrid.u      ,'linewidth',3,'color',c_with,'linestyle','--'); % hybrid
line('parent',H.Ax(2),'xdata',UFEMISM.F0.DIVA_sans.x  ,'ydata',UFEMISM.F0.DIVA_sans.u   ,'linewidth',3,'color',c_sans,'linestyle','-');  % DIVA_sans
line('parent',H.Ax(2),'xdata',UFEMISM.F0.hybrid_sans.x,'ydata',UFEMISM.F0.hybrid_sans.u ,'linewidth',3,'color',c_sans,'linestyle','--'); % hybrid_sans

line('parent',H.Ax(3),'xdata',UFEMISM.F1.DIVA.x       ,'ydata',UFEMISM.F1.DIVA.Hs       ,'linewidth',3,'color',c_with,'linestyle','-');  % DIVA
line('parent',H.Ax(3),'xdata',UFEMISM.F1.hybrid.x     ,'ydata',UFEMISM.F1.hybrid.Hs     ,'linewidth',3,'color',c_with,'linestyle','--'); % hybrid
line('parent',H.Ax(3),'xdata',UFEMISM.F1.DIVA_sans.x  ,'ydata',UFEMISM.F1.DIVA_sans.Hs  ,'linewidth',3,'color',c_sans,'linestyle','-');  % DIVA_sans
line('parent',H.Ax(3),'xdata',UFEMISM.F1.hybrid_sans.x,'ydata',UFEMISM.F1.hybrid_sans.Hs,'linewidth',3,'color',c_sans,'linestyle','--'); % hybrid_sans

line('parent',H.Ax(4),'xdata',UFEMISM.F1.DIVA.x       ,'ydata',UFEMISM.F1.DIVA.u        ,'linewidth',3,'color',c_with,'linestyle','-');  % DIVA
line('parent',H.Ax(4),'xdata',UFEMISM.F1.hybrid.x     ,'ydata',UFEMISM.F1.hybrid.u      ,'linewidth',3,'color',c_with,'linestyle','--'); % hybrid
line('parent',H.Ax(4),'xdata',UFEMISM.F1.DIVA_sans.x  ,'ydata',UFEMISM.F1.DIVA_sans.u   ,'linewidth',3,'color',c_sans,'linestyle','-');  % DIVA_sans
line('parent',H.Ax(4),'xdata',UFEMISM.F1.hybrid_sans.x,'ydata',UFEMISM.F1.hybrid_sans.u ,'linewidth',3,'color',c_sans,'linestyle','--'); % hybrid_sans

legend(H.Ax(1),'DIVA','hybrid','DIVA sans','hybrid sans','Full-Stokes','Higher-Order','FS mean','HO mean','location','northeast')