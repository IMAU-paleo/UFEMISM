function H = PlotMeshData(mesh,d,edgecolor)
  if length(d)==mesh.nV
    H = PlotMeshData_a( mesh,d,edgecolor);
  elseif length(d)==mesh.nTri
    H = PlotMeshData_b( mesh,d,edgecolor);
  elseif length(d)==mesh.nAc
    H = PlotMeshData_c( mesh,d,edgecolor);
  else
    error('whaa')
  end
end