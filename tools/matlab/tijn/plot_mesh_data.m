function H = plot_mesh_data( mesh,d,edgecolor)

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
  
  colorbar(H.Ax);
  set(H.Ax,'position',[25 25 aw ah]);
  set(H.Ax,'units','normalized')

  if     (length( d) == mesh.nV)
    H.Patch = plot_mesh_data_a( gca, mesh, d, edgecolor);
  elseif (length( d) == mesh.nTri)
    H.Patch = plot_mesh_data_b( gca, mesh, d, edgecolor);
  elseif (length( d) == mesh.nAc)
    H.Patch = plot_mesh_data_c( gca, mesh, d, edgecolor);
  else
    error('whaa')
  end
  
end