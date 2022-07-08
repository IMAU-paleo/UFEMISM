function data = ufe_PlotMeshFile(folder_name, domain, p_var, n_slice, ...
                                 n_fig, t_pau, p_range, p_diff)

% ufe_PlotMeshFile Plot a UFEMISM output mesh variable
%    ufe_PlotMeshFile(folder_name,domain,p_var) creates a new figure where
%    the evolution of the variable 'p_var' from the domain 'domain' is
%    plotted based on its time dimension from the mesh output files inside
%    the folder 'folder_name' (help_fields_*). If no value is given for
%    'p_var', the mesh will be plotted. The input arguments 'folder_name'
%    and domain are mandatory.
%
%    Example: ufe_PlotMeshFile('results_test','GRL','Hi');
%
%    The output files (and thus the config file used to run the simulation)
%    must contain the variable 'mask' for the specified domain, which is
%    heavily used by this script.
%
%    ufe_PlotMeshFile(...,n_slice) plots the variable 'p_var' only at the
%    specified time slices 'n_slice'. Valid values for 'n_slice' are 'all',
%    'erst', 'last', or a char vector with a single number (e.g. '-21000').
%    This number must be contained within the start and end times of the
%    simulation, but does not need to coincide with a time slice; in this
%    case the closest time slice to the desired time point 'n_slice' will
%    be plotted instead. If no value is given, or if this argument is 0, it
%    defaults to 'all'.
%
%    ufe_PlotMeshFile(...,~,n_fig) plots 'p_var' using a specific figure
%    number 'n_fig', where 'n_fig' is a positive integer. If no value is
%    given, or if this value is 0, it opens a new figure.
%
%    ufe_PlotMeshFile(...,~,~,t_pau) plots the evolution of 'p_var' with a
%    pause-time of 't_pau' before updating to the next time slice, where
%    't_pau' is a positive number. If no value is given, it defaults to a
%    pause of 0 seconds, which although quite fast still allows to see
%    the animation.
%
%    ufe_PlotMeshFile(...,~,~,~,p_range) plots 'p_var' using a specific
%    range 'p_range' for the colormap/colorbar, where 'p_range' is a two-
%    element vector (e.g. [0 4e3]). If no value is given, or if this value
%    is 0, it defaults to the min and max values of 'p_var' across all time
%    slices, unless custom values for common variables were already defined
%    within this script in the "Default settings" section.
%
%    ufe_PlotMeshFile(...,~,~,~,~,p_diff) plots 'p_var' as its difference
%    w.r.t. the initial value in the mesh file, i.e. the values right after
%    that particular mesh was created. Valid values for 'p_diff' are 1
%    (plot diffs) and 0 (do not). If no value is given, 'p_diff' defaults
%    to 0. Note that UFEMISM can output the variables dHi, dHs, and dHb
%    natively, which can be plotted with this script as well and contain
%    the differences of Hi, Hs, and Hb w.r.t. to present-day, respectively.
%
%    data = ufe_PlotMeshFile(...) returns the struct 'data' containing the
%    values of 'p_var' for the plotted time slices, the corresponding time
%    vector, its mask, and mesh data. If more than one mesh output file is
%    present in 'folder_name' for the domain 'domain', and all time slices
%    are requested, the structure 'data' will have a size of 1xN, where N
%    is the number of mesh files in folder 'folder_name'.
%
%    Example:
%
%    hi_grl = ufe_PlotMeshFile('results_test','GRL','Hi','all',2,.5,0,1);
%
%    The example above will plot all time slices of variable Hi for the
%    Greenland domain on figure number 2, with a pause of .5 s between
%    slices, with the default range ([0 3e3]), and showing the differences
%    w.r.t the first values of its corresponding mesh, and then will store
%    these data in the output variable hi_grl.

% === Input checks ===
% ====================

% Number of input arguments
if nargin < 8
   p_diff = 0;
end
if nargin < 7
   p_range = 0;
end
if nargin < 6
   t_pau = 0.5;
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
if nargin < 2
   error('->[E:] Need at least folder path and desired domain as input!')
end

% Valid values: p_diff
if ~isnumeric(p_diff) || all(p_diff ~= [0,1])
   error(['->[E:] Input argument p_diff must be 0 (no diffs) or '          ...
          '1 (diffs)! Attempted value: ' num2str(p_diff)]);
end
% Valid values: p_range
if p_range ~= 0
   if length(p_range)~=2 || ~isnumeric(p_range) || p_range(2) < p_range(1)
      error(['->[E:] Input argument p_range must be a 1x2 numeric vector'  ...
                   ' in which the second element is greater than the first'...
                   ' element! Attempted value: ' num2str(p_range)]);
   end
end
% Valid values: t_pau
if length(t_pau)~=1 || ~isnumeric(t_pau) || t_pau<0
   error(['->[E:] Input argument t_pau must be a scalar non-negative real' ...
                ' number! Attempted value: ' num2str(t_pau)]);
end
% Valid values: n_fig
if length(n_fig)~=1 || ~isnumeric(n_fig) || n_fig<0
   error(['->[E:] Input argument n_fig must be a scalar positive integer!' ...
                ' Attempted value: ' num2str(n_fig)]);
end
% Valid values: n_slice
if n_slice == 0
   n_slice = 'all';
   disp('->[W]: Plotting all time slices... For time=0, use "0" (string).')
elseif ~any(strcmp(n_slice,{'all','erst','last'}))
   if ~ischar(n_slice)
      error('->[E]: Time slice must be a char type or 0 (== "all")!')
   elseif isnan(str2double(n_slice))
      error('->[E]: Value for input argument n_slice not recognised!')
   end
end

% === Default settings ===
% ========================

   switch p_var
       case {'Hi','Hs'}                % Ice thickness or surface elevation
           if p_range == 0             % If no value range was specified
              if p_diff == 0           % If showing normal values
                 p_range = [0 3e3];    % Use default range
              else                     % If showing differences
                 p_range = [-1e3 1e3]; % Use default range for differences
              end
           end

           if p_diff == 0              % If showing normal values
              cmap = (turbo(256));    % Use default colormap
              mask_ice = 1;            % Only show where mask has ice
              p_trans  = 0;            % Only show where data are not 0
           else                        % If showing differences
              cmap = turbo(20);        % Use custom colormap
              mask_ice = 0;            % Show everywhere
              p_trans  = 0;            % Only show where data are not 0
           end

       case {'Hb'}                     % Bedrock elevation
           if p_range == 0             % If no value range was specified
              p_range = [-3e3 3e3];    % Use default range
           end
           mask_ice = 0;               % Show everywhere
           p_trans  = 0;               % Only show where data are not 0
           cmap = (turbo(256));       % Use default colormap

       case {'SL'}                     % Sea level w.r.t present-day
           if p_range == 0             % If no value range was specified
              p_range = [-120 20];     % Use default range
           end
           mask_ice = -3;              % Only show where mask=ocean/shelves
           p_trans  = Inf;             % No value-based transparency
           cmap = (turbo(256));       % Use default colormap

       case {'dHb'}                    % Bedrock deformation w.r.t. PD
           if p_range == 0             % If no value range was specified
              p_range = [-5e2 5e2];    % Use default range
           end
           mask_ice = 0;               % Show everywhere
           p_trans  = 0;               % Only show where data are not 0
           cmap = (turbo(256));       % Use default colormap

       case {'phi_fric'}               % Till angle for sliding law
           if p_range == 0             % If no value range was specified
             %p_range = [-2.5 47.5];
              p_range = [-1 1.5];      % Use default range
           end
           mask_ice = 2;               % Only show where ice is grounded
           p_trans  = inf;             % No value-based transparency
          %cmap = (flipud(jet(10)));
           cmap = (flipud(jet(256)));  % Use custom colormap

       case {'BMB','BMB_shelf'}        % Basal mass balance
           if p_range == 0             % If no value range was specified
              p_range = [-5 5];        % Use default range
           end
           mask_ice = 3;               % Only show where ice is floating
           p_trans  = Inf;             % No value-based transparency
           cmap = (flipud(jet(256)));  % Use custom colormap

       case {'C_abl_constant_inv', ... % Threshold ablation factor
             'C_abl_Ts_inv' ...        % Temperature ablation factor
             'C_abl_Q_inv', ...        % Insolation ablation factor
             'C_refr_inv'}             % Refreezing factor
           mask_ice = 2;               % Only show where ice is grounded
           p_trans  = 0;               % Only show where data are not 0
           cmap = (turbo(256));       % Use default colormap

       case {'uabs_surf','uabs_base'}  % Ice velocities
           if p_range == 0             % If no value range was specified
              p_range = [-1 3.5];      % Use default (log10) range
           end
           mask_ice = 1;               % Only show where mask has ice
           p_trans  = 0;               % Only show where data are not 0
           cmap = (turbo(256));        % Use custom colormap

       otherwise                       % All other variables
           if p_range == 0             % If no value range was specified
              p_range = [0 0];         % Let the min and max set the range
           end
           mask_ice = 0;               % Show everywhere
           p_trans  = 0;               % Only show where data are not 0
           cmap = (turbo(256));       % Use default colormap
   end

% === Plot meshes/masks/data ===
% ==============================

% If specified folder exists
if exist(folder_name,'dir')
   % Get list of help_fields mesh files
   file_list = dir(fullfile(folder_name,['help_fields_' domain '_*.nc']));
else
   error(['->[E:] Could not find folder "' folder_name '"!'])
end

% Check that we actually have valid files
if isempty(file_list)
   error(['->[E:] No mesh files for domain "' domain '" found in folder "' ...
                  folder_name '"!']);
end

% Time points of interest
if strcmp(n_slice,'last')          % If interested in end of run
   ii_start = length(file_list);   % Only loop over last file on list
   ii_end = length(file_list);     % Only loop over last file on list
   t_slice = NaN;                  % Do not look for specific time slice

elseif strcmp(n_slice,'erst')      % If interested in start of run only
   ii_start = 1;                   % Only loop over first file on list
   ii_end = 1;                     % Only loop over first file on list
   t_slice = NaN;                  % Do not look for specific time slice

elseif strcmp(n_slice,'all')       % If interested in whole run
   ii_start = 1;                   % Loop over all files on list
   ii_end = length(file_list);     % Loop over all files on list
   t_slice = NaN;                  % Do not look for specific time slice

elseif ~isnan(str2double(n_slice)) % If intereted in specific time slice
   ii_start = 1;                   % Loop over all files to look for slice
   ii_end = length(file_list);     % Loop over all files to look for slice
   t_slice = str2double(n_slice);  % Mark time slice of interest
end

% - The big plot loop -
% =====================

% Initialise mesh counter for output argument "data"
kk = 1;

% Loop over all mesh files
for ii = ii_start : ii_end

    % Get name of next file
    file_name = fullfile(folder_name,file_list(ii).name);

    % Read the mesh
    mesh.V    = ncread(file_name,'V');    % Vertices
    mesh.nV   = size(mesh.V,1);           % Number of vertices
    mesh.Tri  = ncread(file_name,'Tri');  % Faces
    mesh.nTri = size(mesh.Tri,1);         % Number of faces
    mesh.R    = ncread(file_name,'R');    % Mesh resolutions

    % Read the time vector
    time = ncread(file_name,'time');

    % Read and format the mask
    try
      mask = ncread(file_name,'mask');    % Read mask
    catch
      error('->[E:] Variable mask not in output file. Ouch!')
    end
    mask(mask==5) = 0;                    % Merge coast and land
    mask(mask==8) = 6;                    % Merge calving front and margin
    mask(mask==2) = 4;                    % Merge lakes and shelves
    mask(mask==3) = 2;                    % Adjust indices (sheet)
    mask(mask==4) = 3;                    % Adjust indices (shelf)
    mask(mask==6) = 4;                    % Adjust indices (margin)
    mask(mask==7) = 5;                    % Adjust indices (grounding line)

    % Get min and max mesh resolution in km, rounded to 1st decimal
    mesh.R_min = round(10*min(mesh.R)/1000)/10;
    mesh.R_max = round(10*max(mesh.R)/1000)/10;

    % - Set variable to be plotted -
    % ==============================

    % Read data
    if strcmp(p_var,'mesh')               % If we want just the mesh
       my_var = zeros(mesh.nV,1);         % Fill in dummy values
       edgecolor = 'k';                   % Edges will be black
       facecolor = 'none';                % Do not color the triangles

    elseif strcmp(p_var,'mask')           % If we want the mask
       my_var = mask;                     % Just use the mask
       edgecolor = 'flat';                % Edge colors follow the mask
       facecolor = 'none';                % Do not color the triangles

    else                                  % If we want any other variable
       try
         my_var = ncread(file_name,p_var);  % Read desired variable
       catch
         error(['->[E:] Variable "' p_var '" not in output file. Ouch!'])
       end
       info = ncinfo(file_name,p_var);    % Get variable info
       edgecolor = 'interp';              % Edges follow interpolated data
       facecolor = 'none';                % Triangles will not be filled
       i_data = my_var(:,1);              % Save first value for reference
    end

    % Time slice of interest (if any)
    t_found = 0;
    if ~isnan(t_slice)
       if t_slice < time(1) || t_slice > time(end)
          if ii == ii_end
             error('->[E]: Slice time outside of simulation range!')
          else
             continue % Specified time not in this mesh file. Try next one.
          end
       else
          % Get slice closest to time of interest
          t_diff = time - t_slice;
          slice = find(abs(t_diff)==min(abs(t_diff)));
          if min(abs(t_diff)) > 0
             disp('->[W]: Desired time not found. Using closest slice.')
          end
          t_found = 1;
       end
    end

    % Time points of interest
    if strcmp(n_slice,'last')       % If interested in end of run only
       my_var = my_var(:,end);      % Keep only last time slice
       if strcmp(p_var,'mesh')      % If plotting the mesh only
          time = time(1);           % Get time of mesh creation
       else                         % If plotting any other variable
          time = time(end);         % Set time accordingly
       end
       mask = mask(:,end);          % Set mask accordingly

    elseif strcmp(n_slice,'erst')   % If interested in start of run only
       my_var = my_var(:,1);        % Keep only first time slice
       time = time(1);              % Set time accordingly
       mask = mask(:,1);            % Set mask accordingly

    elseif strcmp(n_slice,'all')    % If interested in whole run
       my_var = my_var(:,:);        % Keep all time slices
       time = time(:);              % Set time accordingly
       mask = mask(:,:);            % Set mask accordingly

    elseif ~isnan(t_slice)          % if interested in specific time
       my_var = my_var(:,slice);    % Keep only slice closest to time
       if strcmp(p_var,'mesh')      % If plotting the mesh only
          time = time(1);           % Get time of mesh creation
       else                         % If plotting any other variable
          time = time(slice);       % Set time accordingly
       end
       mask = mask(:,slice);        % Set mask accordingly
    end

    % - Initialise figure -
    % =====================

    if n_fig > 0                    % If figure number was specified
       H.Fig = figure(n_fig);       % Use desired figure number
       clf;                         % Clean previous axes
    else                            % If no figure nr. was specified
       if ii == ii_start            % If first iteration
          H.Fig = figure;           % Open new figure
       else                         % If not first iteration
          clf;                      % Simply clean previous axes
       end
    end
    if ii == ii_start               % If a new figure was created
       pause(.1);                   % Give it time to format fig on screen
    end

    % Initialise the axes
    H.Ax  = axes('xtick',[],'ytick',[], ...   % No ticks please
                 'box','on');                 % Use box mode

    % - Plot the evolution of the variable -
    % ======================================

    % Loop over all time slices within this mesh file
    for jj = 1: size(my_var,2)

        % Compute the log10 of data for velocities or friction
        if any(strcmp(p_var,{'uabs_surf','uabs_base','phi_fric'}))
           p_data = log10(my_var(:,jj));
        else
           p_data = my_var(:,jj);
        end

        % If so desired, get diffs of var w.r.t. first values for thi* mesh
        if p_diff == 1 && ~any(strcmp(p_var,{'mask','mesh'}))
           p_data = p_data - i_data;
        end

        % - Patch plot -
        % ==============

        % Plot the mesh in a grey tone, to serve as a base for the plot
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
                        'LineWidth',1.0, ...
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
          cb = colorbar(H.Ax,'FontSize',16);
          set(cb,'YTick',[0:5],'YTickLabel',{'Land','Ocean','Ice sheet', ...
                                             'Ice shelf','Ice margin', ...
                                             'Grounding line','Calving front'});

          % Title
          p_title = ['Model mask (at year: ' num2str(time(jj)) ')'];

        % The mesh
        % ========

        elseif strcmp(p_var,'mesh')
          % Title
          p_title = ['Model mesh ( last update at year: '                  ...
                     num2str(time(jj)) ' )'];

        % All other variables
        % ===================

        else
          % Title
          if p_diff == 0
             p_title = [info.Attributes(1).Value ' ['                      ...
                        info.Attributes(2).Value '] at year: '             ...
                        num2str(time(jj))];
          else
             p_title = [info.Attributes(1).Value ' difference ['           ...
                        info.Attributes(2).Value '] at year: '             ...
                        num2str(time(jj))];
          end
          if any(strcmp(p_var,{'uabs_surf','uabs_base'}))
             p_title = [p_title '  (log10)'];                              %#ok<AGROW>
          end

          % Colormap
          colormap(cmap) % Use default/custom colormap for this var

          % Colorbar
          cb = colorbar(H.Ax,'FontSize',16);

          % User-defined a range
          if any(p_range ~= 0)
             caxis(p_range);
          % Undefined range, not for differences
          elseif p_diff == 0
             if min(my_var(:)) ~= max(my_var(:))
                % Use overall min and max across all time slices
                caxis([min(my_var(:)) max(my_var(:))]);
             end
          % Undefined range, for differences
          elseif p_diff == 1
             if min(p_data(:)) ~= max(p_data(:))
                % Use min and max of current differences
                caxis([min(p_data(:)) max(p_data(:))]);
             end
          % else just use the range already defined by the plot function
          end

          % Initialise transparency
          FaceVertexAlphaData = ones(size(p_data));

          % Mask-based transparency (*excludes margins*)
          if mask_ice == 3
             FaceVertexAlphaData(mask(:,jj)~=3) = 0; % Shelves
          elseif mask_ice == 2
             FaceVertexAlphaData(mask(:,jj)~=2) = 0; % Sheet
             FaceVertexAlphaData(mask(:,jj)==5) = 1; % GL
          elseif mask_ice == 1
             FaceVertexAlphaData(mask(:,jj)<=1) = 0; % Ice
          elseif mask_ice == -1
             FaceVertexAlphaData(mask(:,jj)>=2) = 0; % Ocean/Land
          elseif mask_ice == -2
             FaceVertexAlphaData(mask(:,jj)~=1) = 0; % Ocean
          elseif mask_ice == -3
             FaceVertexAlphaData(mask(:,jj)~=1) = 0; % Ocean
             FaceVertexAlphaData(mask(:,jj)==3) = 1; % Shelves
          end

          % Value-based transparency
          FaceVertexAlphaData(abs(p_data-p_trans)<1e-6) = 0;

## Not yet implemented in Octave :( ############################################
##
##          % Apply transparency
##          set(H.Patch,'FaceVertexAlphaData',FaceVertexAlphaData);
##          set(H.Patch,'EdgeAlpha','interp');
##          set(H.Patch,'FaceAlpha','interp');
##
################################################################################

        end

        % Title
        p_title = [p_title newline 'Resolution: ' num2str(mesh.R_min)      ...
                   ' - ' num2str(mesh.R_max) ' km'];
        tt = title(p_title,'FontSize',18);

        % Axes box
        axis equal tight  % Make axes proportional and tight around data

        % Pause time
        pause(t_pau)      % Pause it a bit before updating to next mesh

        % Delete plotted data so the next iteration has a
        % clean plate (except for final slice, of course)
        if jj < size(my_var,2)
           delete(M.Patch)
           delete(H.Patch)
        end

    end % 1: size(my_var,2)

    % Output data for this mesh
    data(kk).time = time;                                                  %#ok<AGROW>
    data(kk).mask = mask;                                                  %#ok<AGROW>
    data(kk).vals = my_var;                                                %#ok<AGROW>
    data(kk).mesh = mesh;                                                  %#ok<AGROW>

    % Update mesh plot counter
    kk = kk+1;

    % If plotting a specific time slice, exit loop after doing it
    if t_found == 1
       break
    end

end % 1 : length(file_list)
