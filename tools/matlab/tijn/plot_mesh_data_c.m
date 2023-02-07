function H = plot_mesh_data_c( mesh, d_c, edgecolor)

  xr = max(mesh.V(:,1)) - min(mesh.V(:,1));
  yr = max(mesh.V(:,2)) - min(mesh.V(:,2));
  ah = 400;
  aw = ah*xr/yr;
  
  if aw>1200
    ah = ah*1200/aw;
    aw = 1200;
  end
  
  fh = ah+60;
  fw = aw+130;

  H.Fig = figure('position',[900 500 fw fh],'color','w');
  H.Ax  = axes('units','pixels','position',[25 25 aw ah],'xtick',[],'ytick',[],'box','on',...
    'xlim',[min(mesh.V(1:mesh.nV,1)) max(mesh.V(1:mesh.nV,1))],'ylim',[min(mesh.V(1:mesh.nV,2)) max(mesh.V(1:mesh.nV,2))],...
    'fontsize',24);
  H.Patch = patch('parent',H.Ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
    'facevertexcdata',[],'facecolor','none','edgecolor',edgecolor);
  
  colorbar(H.Ax);
  set(H.Ax,'position',[25 25 aw ah]);
  set(H.Ax,'units','normalized')
  
  ncols = 256;
  cmap = parula(ncols);
  clim = [min(d_c), max(d_c)];
  if (clim(1)==clim(2))
    clim(1) = clim(1)-1;
    clim(2) = clim(2)+1;
  end
  set(H.Ax,'clim',clim)
  
  for ci = 1: ncols
    linedata(ci).x = [];
    linedata(ci).y = [];
  end
  for aci = 1: mesh.nAc
    ci = (d_c( aci) - clim(1)) / (clim(2) - clim(1));
    ci = max(1,min(ncols,1 + round(ci*(ncols-1))));
    linedata(ci).x( end+1) = mesh.VAc( aci,1);
    linedata(ci).y( end+1) = mesh.VAc( aci,2);
  end
  for ci = 1: ncols
    line('parent',H.Ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
      'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',8);
  end
  
  % lastly NaN values
  line('parent',H.Ax,'xdata',mesh.VAc( isnan(d_c),1),'ydata',mesh.VAc( isnan(d_c),2),'linestyle','none',...
      'marker','x','markerfacecolor','r','markeredgecolor','r','markersize',8);
  
  drawnow('update');

end