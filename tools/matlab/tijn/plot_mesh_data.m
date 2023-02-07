function H = plot_mesh_data( mesh,d,edgecolor)
  if length(d)==mesh.nV
    H = plot_mesh_data_a( mesh,d,edgecolor);
  elseif length(d)==mesh.nTri
    H = plot_mesh_data_b( mesh,d,edgecolor);
  elseif length(d)==mesh.nAc
    H = plot_mesh_data_c( mesh,d,edgecolor);
  else
    error('whaa')
  end
end