function H = plot_mesh_data_b( ax, mesh, d_b, edgecolor)

H = patch('parent',ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
  'facevertexcdata',d_b,'facecolor','flat','edgecolor',edgecolor);

end