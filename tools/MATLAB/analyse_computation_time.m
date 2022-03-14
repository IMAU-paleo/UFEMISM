function analyse_computation_time( filename)
  
filename = '../../results_20220314_001/resource_tracking.nc';

R = read_resource_tracking_file( filename);

% Set up GUI
wa = 1200;
ha = 800;

wf = 25 + wa + 25;
hf = 25 + ha + 25;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = axes('units','pixels','position',[25,25,wa,ha],'xtick',[],'ytick',[],'fontsize',24);
H.Ax.XAxis.Visible = 'off';
H.Ax.YAxis.Visible = 'off';

% Find maximum rank ever
max_rank = find_max_rank( R);
set(H.Ax,'xlim',[0,max_rank+3]);

% Draw all the patches
R = draw_patches( H, R);

% Hide all except the first rank
for ci = 1: length( R.children)
  hide_all_children( R.children{ ci})
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
      str = [num2str(round(R.children{ ci}.tcomp_rel*1000)/10) ' % - ' R.children{ ci}.name];
      R.children{ ci}.text = text(xt,yt,treat_underscore( str),'fontsize',24,'color',colors(ci,:));
      
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
  function click_on_patch( src,~,tree)
    
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