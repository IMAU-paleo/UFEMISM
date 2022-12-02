function H = plot_mesh_data_c( ax, mesh, d_c, edgecolor)

  H = patch('parent',ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
    'facevertexcdata',[],'facecolor','none','edgecolor',edgecolor);
  
  ncols = 256;
  cmap = parula(ncols);
  clim = [min(d_c), max(d_c)];
  if (clim(1)==clim(2))
    clim(1) = clim(1)-1;
    clim(2) = clim(2)+1;
  end
  set(ax,'clim',clim)
  
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
    line('parent',ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
      'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',8);
  end
  
  % lastly NaN values
  line('parent',ax,'xdata',mesh.VAc( isnan(d_c),1),'ydata',mesh.VAc( isnan(d_c),2),'linestyle','none',...
      'marker','x','markerfacecolor','r','markeredgecolor','r','markersize',8);

end