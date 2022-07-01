function data = ufe_PlotMeshFile(folder_name, domain, p_var, n_slice, n_fig, p_range, t_pau, p_diff)

% === Input checks ===
% ====================

if nargin < 8
   p_diff = 0;
end
if nargin < 7
   t_pau = .05;
end
if nargin < 6
   p_range = 0;
end
if nargin < 5
   n_fig = 0;
end
if nargin < 4
   n_slice = 'all';
end
if nargin < 3
   p_var = 'mesh';
end

% === Default settings ===
% ========================

   switch p_var
       case {'Hi','Hs'}          % Ice thickness or surface elevation
           if p_range == 0
              if p_diff == 0
                 p_range = [0 3e3];
              else
                 p_range = [-1e3 1e3];
              end
           end
           mask_ice = 1;         % Only show where there is ice
           p_trans  = 0;         % Only show where data are not 0

           if p_diff == 0
              cmap = (jet(256));
           else
              cmap = turbo(20);
           end

       case {'Hb'}               % Bedrock elevation
           if p_range == 0
              p_range = [-3e3 3e3];
           end
           mask_ice = 0;         % Show everywhere
           p_trans  = 0;         % Only show where data are not 0
           cmap = (jet(256));

       case {'SL'}               % Sea level w.r.t present-day
           if p_range == 0
              p_range = [-120 20];
           end
           mask_ice = -3;        % Only show where there is ocean/shelves
           p_trans  = Inf;       % No value-based transparency
           cmap = (jet(256));

       case {'dHb'}              % bedrock deformation w.r.t. present-day
           if p_range == 0
              p_range = [-5e2 5e2];
           end
           mask_ice = 0;         % Show everywhere
           p_trans  = 0;         % Only show where data are not 0
           cmap = (jet(256));

       case {'phi_fric'}         % Till angle for sliding law
           if p_range == 0
              p_range = [-2.5 47.5];
           end
           mask_ice = 2;         % Only show where ice is grounded
           p_trans  = 0;         % Only show where data are not 0
           cmap = (flipud(jet(10)));

       case {'BMB','BMB_shelf'}  % Basal mass balance
           if p_range == 0
              p_range = [-5 5];
           end
           mask_ice = 3;         % Only show where there is ocean/shelves
           p_trans  = Inf;       % No value-based transparency
           cmap = (flipud(jet(256)));

       case {'C_abl_constant_inv','C_abl_Ts_inv'}  % IMAU-ITM parameters
           mask_ice = 2;         % Only show where ice is grounded
           p_trans  = 0;         % Only show where data are not 0
           cmap = (jet(256));

       case {'uabs_surf','uabs_base'}  % Ice velocities
           if p_range == 0
              p_range = [-1 3.5];
           end
           mask_ice = 1;         % Only show where there is ice
           p_trans  = 0;         % Only show where data are not 0
           cmap = (turbo(256));

       otherwise                 % All other variables
           if p_range == 0
              p_range = [0 0];
           end
           mask_ice = 1;         % Only show where there is ice
           p_trans  = 0.0227;    % Only show where data are not 0
           cmap = (jet(256));
   end

% === Plot meshes/masks/data ===
% ==============================

% Get list of files
file_list = dir(fullfile(folder_name,['help_fields_' domain '_*.nc']));

% Time points of interest
if strcmp(n_slice,'last')        % If interested in end of run
   ii_start = length(file_list); % Only loop over last file on list
   ii_end = length(file_list);
elseif strcmp(n_slice,'erst')    % If interested in start of run only
   ii_start = 1;                 % Only loop over first file on list
   ii_end = 1;
else                             % If interested in whole run
   ii_start = 1;                 % Loop over all files on list
   ii_end = length(file_list);
end

% - Loop over all mesh files -
% ============================

for ii = ii_start : ii_end

    % Get name of next file
    file_name = fullfile(folder_name,file_list(ii).name);

    % Read mesh and time
    mesh = ReadMeshFromFile(file_name);
    time = ncread(file_name,'time');

    % The mask
    mask = ncread(file_name,'mask');      % Read mask
    mask(mask==5) = 0;                    % Merge coast and land into one
    mask(mask==8) = 6;                    % Merge coast and land into one
    mask(mask==2) = 4;                    % Merge lakes and shelves into one
    mask(mask==3) = 2;                    % Move the indices down (sheet)
    mask(mask==4) = 3;                    % Move the indices down (shelf)
    mask(mask==6) = 4;                    % Move the indices down (margin)
    mask(mask==7) = 5;                    % Move the indices down (grounding line)

    % Get min and max mesh resolution in km, rounded to 1st decimal
    mesh.R_min = round(10*min(mesh.R)/1000)/10;
    mesh.R_max = round(10*max(mesh.R)/1000)/10;

    % - Set variable to be plotted -
    % ==============================

    if strcmp(p_var,'mesh')               % If we want just the mesh
       my_var = zeros(mesh.nV,1);         % Fill in dummy values
       edgecolor = 'k';                   % Edges will be black
       facecolor = 'none';                % Do not color the triangles

    elseif strcmp(p_var,'mask')           % If we want the mask
       my_var = mask;                     % Just use the mask
       edgecolor = 'flat';                % Edge colors will follow the mask
       facecolor = 'none';                % Do not color the triangles

    else                                  % If we want any other variable
       my_var = ncread(file_name,p_var);  % Read desired variable
       info = ncinfo(file_name,p_var);    % Get variable info (name, units, etc.)
       edgecolor = 'interp';              % Edges will follow interpolated data
       facecolor = 'none';                % Triangles will not be filled
       i_data = my_var(:,1);              % Save initial value for reference
    end

    % Time points of interest
    if strcmp(n_slice,'last')             % If interested in end of run only
       my_var = my_var(:,end);            % Keep only last time slice
       if strcmp(p_var,'mesh')
          time = time(1);                 % Get time of mesh creation
       else
          time = time(end);               % Set time accordingly
       end
       mask = mask(:,end);
    elseif strcmp(n_slice,'erst')         % If interested in start of run only
       my_var = my_var(:,1);              % Keep only first time slice
       time = time(1);                    % Set time accordingly
       mask = mask(:,1);
    end

    % - Initialise figure -
    % =====================

    if n_fig > 0                          % If figure number was specified
       H.Fig = figure(n_fig);             % Use desired figure number
       clf;                               % Clean previous axes
    else                                  % If no figure nr. was specified
       if ii == ii_start                  % If first iteration
          H.Fig = figure;                 % Open new figure
       else                               % If not first iteration
          clf;                            % Simply clean previous axes
       end
    end
    if ii == ii_start                     % If a new figure was created
       pause(.1);   % Give it some time to properly format it on screen
    end

    % Initialise the axes
    H.Ax  = axes('xtick',[],'ytick',[], ...
                 'box','on', ...
                 'fontsize',12);

    % - Plot the evolution of the variable -
    % ======================================

    for jj = 1: size(my_var,2)

        % Compute the log10 of the data when plotting velocities
        if any(strcmp(p_var,{'uabs_surf','uabs_base'}))
           p_data = log10(my_var(:,jj));
        else
           p_data = my_var(:,jj);
        end

        % Get difference of variable w.r.t. initial values on *this* mesh
        if p_diff == 1 && ~any(strcmp(p_var,{'mask','mesh'}))
           p_data = p_data - i_data;
        end

        % - Patch plot -
        % ==============

        % Plot first the mesh in a grey tone, to serve as a base for the plot
        M.Patch = patch('Parent',H.Ax, ...
                        'Vertices',mesh.V(1:mesh.nV,:), ...
                        'Faces', mesh.Tri(1:mesh.nTri,:), ...
                        'FaceVertexCData',zeros(mesh.nV,1), ...
                        'FaceColor',facecolor, ...
                        'EdgeColor',[.8 .8 .8], ...
                        'LineWidth',1.0, ...
                        'LineStyle','-');

        % Keep the mesh in place; next plot will go on top of it
        hold on

        % Plot the mesh and (optionally) another quantity
        H.Patch = patch('Parent',H.Ax, ...
                        'Vertices',mesh.V(1:mesh.nV,:), ...
                        'Faces', mesh.Tri(1:mesh.nTri,:), ...
                        'FaceVertexCData',p_data, ...
                        'FaceColor',facecolor, ...
                        'EdgeColor',edgecolor, ...
                        'LineWidth',1.5, ...
                        'LineStyle','-');

        % - Cosmetics -
        % =============

        % The mask
        % ========

        if strcmp(p_var,'mask')
          % Colormap
          cmap = [[211, 182, 131]/255; ... % Land
                  [  3, 113, 156]/255; ... % Ocean
                  [152, 239, 249]/255; ... % Ice sheet
                  [ 60, 153, 146]/255; ... % Ice shelf
                  [  0,   0,   0]/255; ... % Ice Margin
                  [254,  44,  84]/255];    % Grounding line
          colormap(cmap);
          caxis([-.5 5.5]);
          % Colorbar
          cb = colorbar(H.Ax);
          cb.TickLabels = {'Land','Ocean','Ice sheet','Ice shelf', ...
                           'Ice margin','Grounding line','Calving front'};
          cb.FontSize = 16;
          % Title
          p_title = ['Model mask (at year: ' num2str(time(jj)) ')'];

        % The mesh
        % ========

        elseif strcmp(p_var,'mesh')
          % Title
          p_title = ['Model mesh ( last update at year: ' num2str(time(jj)) ' )'];

        % All other variables
        % ===================

        else
          % Title
          if p_diff == 0
             p_title = [info.Attributes(1).Value ' [' ...
                        info.Attributes(2).Value '] at year: ' ...
                        num2str(time(jj))];
          else
             p_title = [info.Attributes(1).Value ' difference [' ...
                        info.Attributes(2).Value '] at year: ' ...
                        num2str(time(jj))];
          end
          if any(strcmp(p_var,{'uabs_surf','uabs_base'}))
             p_title = [p_title '  (log10)'];                              %#ok<AGROW>
          end
          % Colormap
          colormap(cmap)
          % Colorbar
          cb = colorbar(H.Ax);
%         cb.FontSize = 16;
          % Data range
          if any(p_range ~= 0)
             caxis(p_range);
          elseif min(my_var(:)) ~= max(my_var(:))
             caxis([min(my_var(:)) max(my_var(:))]);
          end
          % Value-based transparency
 %         H.Patch.FaceVertexAlphaData = double(abs(my_var(:,jj)-p_trans)>1e-12);
          % Mask-based transparency
          if mask_ice == 3
 %            H.Patch.FaceVertexAlphaData(mask(:,jj)~=3) = 0; % Shelves
            %H.Patch.FaceVertexAlphaData(mask(:,jj)==5) = 1; % GL
          elseif mask_ice == 2
 %            H.Patch.FaceVertexAlphaData(mask(:,jj)~=2) = 0; % Sheet
 %            H.Patch.FaceVertexAlphaData(mask(:,jj)==5) = 1; % GL
          elseif mask_ice == 1
 %            H.Patch.FaceVertexAlphaData(mask(:,jj)<=1) = 0; % Ice
          elseif mask_ice == -1
 %            H.Patch.FaceVertexAlphaData(mask(:,jj)>=2) = 0; % Ocean/Land
          elseif mask_ice == -2
 %            H.Patch.FaceVertexAlphaData(mask(:,jj)~=1) = 0; % Ocean
          elseif mask_ice == -3
 %            H.Patch.FaceVertexAlphaData(mask(:,jj)~=1) = 0; % Ocean
 %            H.Patch.FaceVertexAlphaData(mask(:,jj)==3) = 1; % Shelves
            %H.Patch.FaceVertexAlphaData(mask(:,jj)==4) = 1; % Margins
          end
          % Apply transparency
 %         H.Patch.EdgeAlpha = 'interp';
 %         H.Patch.FaceAlpha = 'interp';
        end

        % Title
        tt = title(p_title);
%        tt.FontSize = 20;
        % Subtitle
%        st = subtitle(['Resolution: ' num2str(mesh.R_min) ...
%                                ' - ' num2str(mesh.R_max) ' km']);
%        st.FontSize = 18;
        % Final touches
        axis equal tight  % Make axes proportional and tight around data
        pause(t_pau)      % Pause it a bit before updating to next mesh

        % Delete plotted data so the next iteration has a
        % clean plate (except for final slice, of course)
        if jj < size(my_var,2)
           delete(M.Patch)
           delete(H.Patch)
        end

    end % 1: size(my_var,2)

end % 1 : length(file_list)l

% Output
data.time = time;
data.mask = mask;
data.vals = my_var;
data.mesh = mesh;
