function H = plot_mesh_data_a( ax, mesh, d_a, edgecolor)

H = patch('parent',ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
  'facevertexcdata',d_a,'facecolor','none','edgecolor','interp','linewidth',3);
  
end