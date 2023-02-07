function H = plot_mesh_data_b( mesh, d_b, edgecolor)

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
    'facevertexcdata',d_b,'facecolor','flat','edgecolor',edgecolor);
  
  colorbar(H.Ax);
  set(H.Ax,'position',[25 25 aw ah]);
  set(H.Ax,'units','normalized')
  
  drawnow('update');

end