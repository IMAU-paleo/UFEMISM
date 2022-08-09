function analyse_computation_time( filename)
  
filename = '../../results_20220315_001/resource_tracking.nc';

% Set up GUI
wa = 1000;
ha = 800;
wa2 = 400;
ha2 = ((ha+50)-100-25-25)/2;

wf = 25 + wa + 25 + wa2 + 100;
hf = 25 + ha + 25;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = axes('units','pixels','position',[25,25,wa,ha],'xtick',[],'ytick',[],'fontsize',24);
H.Ax.XAxis.Visible = 'off';
H.Ax.YAxis.Visible = 'off';

x = 25 + wa + 25;
y1 = 100;
y2 = 100 + ha2 + 25;
H.Ax2 = axes('units','pixels','position',[x,y2,wa2,ha2],'fontsize',24,'xgrid','on','ygrid','on','yaxislocation','right','xticklabels',[]);
H.Ax3 = axes('units','pixels','position',[x,y1,wa2,ha2],'fontsize',24,'xgrid','on','ygrid','on','yaxislocation','right');
xlabel(H.Ax3,'Model time (yr)')
ylabel(H.Ax2,'Computation time (s)')
ylabel(H.Ax3,'Memory use (MB)')
H.line_tcomp = line('parent',H.Ax2,'xdata',[],'ydata',[],'linewidth',3);
H.line_mem   = line('parent',H.Ax3,'xdata',[],'ydata',[],'linewidth',3);

% Read resource tracking file
R = read_resource_tracking_file( filename);

% Find maximum rank ever
max_rank = find_max_rank( R);
set(H.Ax,'xlim',[0,max_rank+3]);

% Draw all the patches
R = draw_patches( H, R);

% Hide all except the first rank
for ci = 1: length( R.children)
  hide_all_children( R.children{ ci})
end

  function R = read_resource_tracking_file( filename)
    % Read a resource tracking file generated by IMAU-ICE or UFEMISM

    %% Read data from NetCDF file
    R = [];
    nr = 0;

    for vi = 1:1000

      if vi < 10
        vi_str = ['0000' num2str(vi)];
      elseif vi < 100
        vi_str = ['000' num2str(vi)];
      elseif vi < 1000
        vi_str = ['00' num2str(vi)];
      elseif vi < 10000
        vi_str = ['0' num2str(vi)];
      elseif vi < 100000
        vi_str = num2str(vi);
      end

      var_name_name  = ['name_'  vi_str];
      var_name_tcomp = ['tcomp_' vi_str];
      var_name_mem   = ['mem_'   vi_str];

      % Read subroutine name
      path_int_enc = ncread( filename, var_name_name);
      routine_path = decode_path_int( path_int_enc);
      if strcmp(routine_path( 1:length('subroutine_placeholder')),'subroutine_placeholder'); break; end

      % Read model time
      time  = ncread( filename, 'time');

      % Read computation time
      tcomp = ncread( filename, var_name_tcomp);

      % Read memory use
      mem   = ncread( filename, var_name_mem);

      % Add to data structure
      nr = nr+1;
      R( nr).routine_path = routine_path;
      R( nr).time         = time;
      R( nr).tcomp        = tcomp;
      R( nr).tcomp_tot    = sum( tcomp);
      R( nr).mem          = mem;

    end

    %% Organise the resource use hierarchically
    R = organise_resource_hierarchy( R);

    %% Add the root parent
    tcomp_tot = 0;
    for ri = 1:length(R)
      tcomp_tot = tcomp_tot + R( ri).tcomp_tot;
    end

    R_old = R;
    R     = [];

    R.rank         = 1;
    R.routine_path = 'IMAU_ICE_program';
    R.name         = 'IMAU_ICE_program';
    R.tree{1}      = 'IMAU_ICE_program';
    R.tcomp_tot    = tcomp_tot;
    R.children = {};
    for ci = 1: length(R_old)
      R.children{end+1} = R_old( ci);
    end

    %% Sort resources use
    R = sort_resources( R);
    R.mem = 0;
    R = sum_childrens_memory( R);
    R = calc_relative_total_computation_time( R, R.tcomp_tot);

    function routine_path = decode_path_int( path_int_enc)
      % Convert the routine path from integr form back to a string

      routine_path = blanks(1024);
      for i = 1:length( path_int_enc)
        switch path_int_enc( i)
          case -1
            break
          case 0
            routine_path( i) = '0';
          case 1
            routine_path( i) = '1';
          case 2
            routine_path( i) = '2';
          case 3
            routine_path( i) = '3';
          case 4
            routine_path( i) = '4';
          case 5
            routine_path( i) = '5';
          case 6
            routine_path( i) = '6';
          case 7
            routine_path( i) = '7';
          case 8
            routine_path( i) = '8';
          case 9
            routine_path( i) = '9';
          case 11
            routine_path( i) = 'a';
          case 12
            routine_path( i) = 'b';
          case 13
            routine_path( i) = 'c';
          case 14
            routine_path( i) = 'd';
          case 15
            routine_path( i) = 'e';
          case 16
            routine_path( i) = 'f';
          case 17
            routine_path( i) = 'g';
          case 18
            routine_path( i) = 'h';
          case 19
            routine_path( i) = 'i';
          case 20
            routine_path( i) = 'j';
          case 21
            routine_path( i) = 'k';
          case 22
            routine_path( i) = 'l';
          case 23
            routine_path( i) = 'm';
          case 24
            routine_path( i) = 'n';
          case 25
            routine_path( i) = 'o';
          case 26
            routine_path( i) = 'p';
          case 27
            routine_path( i) = 'q';
          case 28
            routine_path( i) = 'r';
          case 29
            routine_path( i) = 's';
          case 30
            routine_path( i) = 't';
          case 31
            routine_path( i) = 'u';
          case 32
            routine_path( i) = 'v';
          case 33
            routine_path( i) = 'w';
          case 34
            routine_path( i) = 'x';
          case 35
            routine_path( i) = 'y';
          case 36
            routine_path( i) = 'z';
          case 37
            routine_path( i) = 'A';
          case 38
            routine_path( i) = 'B';
          case 39
            routine_path( i) = 'C';
          case 40
            routine_path( i) = 'D';
          case 41
            routine_path( i) = 'E';
          case 42
            routine_path( i) = 'F';
          case 43
            routine_path( i) = 'G';
          case 44
            routine_path( i) = 'H';
          case 45
            routine_path( i) = 'I';
          case 46
            routine_path( i) = 'J';
          case 47
            routine_path( i) = 'K';
          case 48
            routine_path( i) = 'L';
          case 49
            routine_path( i) = 'M';
          case 50
            routine_path( i) = 'N';
          case 51
            routine_path( i) = 'O';
          case 52
            routine_path( i) = 'P';
          case 53
            routine_path( i) = 'Q';
          case 54
            routine_path( i) = 'R';
          case 55
            routine_path( i) = 'S';
          case 56
            routine_path( i) = 'T';
          case 57
            routine_path( i) = 'U';
          case 58
            routine_path( i) = 'V';
          case 59
            routine_path( i) = 'W';
          case 60
            routine_path( i) = 'X';
          case 61
            routine_path( i) = 'Y';
          case 62
            routine_path( i) = 'Z';
          case 63
            routine_path( i) = '_';
          case 64
            routine_path( i) = '/';
          case 65
            routine_path( i) = '(';
          case 66
            routine_path( i) = ')';
          otherwise
            error('whaa!')
        end
      end
      routine_path = routine_path( 1:i-1);
    end
    function R = organise_resource_hierarchy( R)
      % Organise the resource use hierarchically

      % Get the rank and tree for each subroutine
      max_rank = 0;
      for ri = 1: length( R)

        rank = 0;
        tree = [];

        pathstr = R( ri).routine_path;

        while ~isempty( pathstr)
          ci = 1;
          while ci <= length( pathstr)
            if strcmp( pathstr( ci),'/'); break; end
            ci = ci+1;
          end
          rank = rank+1;
          tree{ end+1} = pathstr( 1:ci-1);
          pathstr = pathstr( ci+1:end);
        end

        max_rank = max( max_rank, rank);

        R( ri).name = tree{ end};
        R( ri).rank = rank;
        R( ri).tree = tree;
        R( ri).children = {};

      end


      % For each subroutine, see if we can find its parent. If so, add it there.
      for rank = max_rank:-1:1

        ri = 1;
        while ri <= length( R)

          if R( ri).rank ~= rank
            ri = ri + 1;
            continue
          end

          found_parent = false;
          for rj = 1: length( R)
            if is_parent( R( ri), R( rj))
              found_parent = true;
              break
            end
          end
          if found_parent
            R( rj).children{ end+1} = R( ri);
            R( ri) = [];
          else
            ri = ri+1;
          end

        end

      end

      function isso = is_parent( sub_child, sub_parent)
        isso = false;
        if sub_parent.rank ~= sub_child.rank - 1; return; end
        for j = 1: sub_parent.rank
          if ~strcmp( sub_parent.tree{ j}, sub_child.tree{ j}); return; end
        end
        isso = true;
      end

    end
    function R = sort_resources( R)

      tcomp_tot = zeros( length(R.children),1);
      for ci = 1: length(R.children)
        tcomp_tot( ci) = R.children{ ci}.tcomp_tot;
      end

      [~,ind] = sortrows( tcomp_tot);

      children_old = R.children;
      for ci = 1: length(R.children)
        R.children{ci} = children_old{ ind( ci)};
      end

      for ci = 1: length(R.children)
        R.children{ ci} = sort_resources( R.children{ ci});
      end

    end
    function R = sum_childrens_memory( R)
      % set R.mem to the sum of both R's own memory use and that of all its
      % children (recursive)

      if isempty(R.children)
        return
      else
        mem_children = 0;
        for cci = 1: length( R.children)
          R.children{ cci} = sum_childrens_memory( R.children{ cci});
          mem_children = mem_children + R.children{ cci}.mem;
        end
        R.mem = R.mem + mem_children;
      end
    end
    function R = calc_relative_total_computation_time( R, tcomp_tot)
      % For each subroutine, calculate its computation time relative to the
      % total model time (recursive)

      R.tcomp_tot_rel = R.tcomp_tot / tcomp_tot;
      for cci = 1: length( R.children)
        R.children{ cci} = calc_relative_total_computation_time( R.children{ cci}, tcomp_tot);
      end

    end

  end
  function R = draw_patches( H, R)
    
    rank = R.rank;
    
    % Calculate relative computation time for all children
    for ci = 1: length( R.children)
      R.children{ ci}.tcomp_rel = R.children{ ci}.tcomp_tot / R.tcomp_tot;
    end
    
    % Draw patches and text labels for all children
    w  = 0.8;
    xl = (1-w)/2 + (rank-1);
    xu = xl + w;
    yl = 0;
    
    colors = flipud(lines( length( R.children)));
    
    for ci = 1: length( R.children)
      yu = yl + R.children{ ci}.tcomp_rel;
      V = [xl,yl;xu,yl;xu,yu;xl,yu];
      R.children{ ci}.patch = patch('parent',H.Ax,'vertices',V,'faces',[1,2,3,4],'edgecolor','k','facecolor',colors(ci,:),...
        'ButtonDownFcn',{@click_on_patch,R.children{ ci}.tree});
      
      % Text labels to the right
      xt = xl + 1;
      yt = (-0.5 + ci)/length(R.children);
      str = [num2str(round(R.children{ ci}.tcomp_rel*1000)/10) ' % ( ' num2str(round(R.children{ ci}.tcomp_tot_rel*1000)/10) '%) - ' R.children{ ci}.name];
      R.children{ ci}.text = text(H.Ax,xt,yt,treat_underscore( str),'fontsize',20,'color',colors(ci,:));
      
      yl = yu;
    end
    
    % Recursive
    for ci = 1: length( R.children)
      R.children{ ci} = draw_patches( H, R.children{ ci});
    end
    
  end
  function max_rank = find_max_rank( resources)
    max_rank = 0;
    for ri = 1: length( resources)
      max_rank = max( max_rank, resources( ri).rank);
      if ~isempty( resources( ri).children)
        for ci = 1: length( resources( ri).children)
          max_rank = max( max_rank, find_max_rank( resources( ri).children{ ci}));
        end
      end
    end
  end
  function str = treat_underscore( str)
    i = 1;
    while i <= length(str)
      if strcmp( str( i),'_')
        str = [str(1:i-1) '\_' str(i+1:end)];
        i = i+2;
      else
        i = i+1;
      end
    end
  end
  function show_direct_children( R)
    for ci = 1: length( R.children)
      set(R.children{ ci}.patch,'visible','on');
      set(R.children{ ci}.text,'visible','on');
    end
  end
  function hide_all_children( R)
    for ci = 1: length( R.children)
      set(R.children{ ci}.patch,'visible','off');
      set(R.children{ ci}.text,'visible','off');
      hide_all_children( R.children{ ci});
    end
  end
  function click_on_patch( ~,~,tree)
    
    % Find the current subroutine in the R structure
    tree(1) = [];
    R_temp = R;
    R_parent = R_temp;
    while ~isempty(tree)
      for ci = 1: length( R_temp.children)
        if strcmpi( R_temp.children{ ci}.name,tree{1})
          R_parent = R_temp;
          R_temp = R_temp.children{ ci};
          tree(1) = [];
          break
        end
      end
    end
    
    % Update the computation time and memory use graphs
    set(H.line_tcomp,'xdata',R_temp.time,'ydata',R_temp.tcomp);
    set(H.line_mem  ,'xdata',R_temp.time,'ydata',R_temp.mem/1e6); % Divide by 1e to show memory use in MB
    
    % Check if this subroutine has children; if not, do nothing
    if isempty(R_temp.children)
      return
    end
    
    % Check if this patch's children are visible. If so, hide them. If not,
    % show them.
    if strcmpi(R_temp.children{1}.patch.Visible,'off')
      % Children are currently hidden; show them
      
      % Hide text for this subroutine and its siblings
      for ci = 1: length( R_parent.children)
        set( R_parent.children{ ci}.text,'visible','off')
      end
      
      % Hide patches and text for this subroutine's siblings' children
      for ci = 1: length( R_parent.children)
        hide_all_children( R_parent.children{ ci})
      end

      % Show this subroutines direct children
      show_direct_children( R_temp)
      
    else
      % Children are currently visible; hide them

      % Hide all this subroutines children
      hide_all_children( R_temp)
      
      % Show text for this subroutine and its siblings
      for ci = 1: length( R_parent.children)
        set( R_parent.children{ ci}.text,'visible','on')
      end
    end
     
  end

end