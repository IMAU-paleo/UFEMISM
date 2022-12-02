function H = plot_mesh( mesh)

  xr = max( mesh.V( :,1)) - min( mesh.V( :,1));
  yr = max( mesh.V( :,2)) - min( mesh.V( :,2));
  ah = 500;
  aw = ah*xr/yr;
  fh = ah+50;
  fw = aw+50;

  H.Fig = figure('position',[900 500 fw fh],'color','w');
  H.Ax  = axes('units','pixels','position',[25 25 aw ah],'xtick',[],'ytick',[],'box','on',...
    'xlim',[min(mesh.V(1:mesh.nV,1)) max(mesh.V(1:mesh.nV,1))],'ylim',[min(mesh.V(1:mesh.nV,2)) max(mesh.V(1:mesh.nV,2))],'fontsize',24);

  H.Patch = patch('parent',H.Ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),'facecolor','none','edgecolor','k');
  drawnow('update');
  
  set(H.Ax,'units','normalized')
end