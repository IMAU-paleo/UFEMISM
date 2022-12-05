function mesh = read_mesh_from_file( filename)

  mesh.V              = ncread(filename,'V');
  mesh.nV             = size(mesh.V,1);
  mesh.Tri            = ncread(filename,'Tri');
  mesh.nTri           = size(mesh.Tri,1);
  
  mesh.xmin           = mesh.V(1,1);
  mesh.xmax           = mesh.V(2,1);
  mesh.ymin           = mesh.V(2,2);
  mesh.ymax           = mesh.V(3,2);
  mesh.tol_dist       = (mesh.xmax - mesh.xmin) * 1e-9;
  
  mesh.nC             = ncread(filename,'nC');
  mesh.C              = ncread(filename,'C');
  mesh.nC_mem         = size( mesh.C,2);
  
  mesh.edge_index     = ncread(filename,'edge_index');
  mesh.Tri_edge_index = ncread(filename,'Tri_edge_index');
  
  mesh.niTri          = ncread(filename,'niTri');
  mesh.iTri           = ncread(filename,'iTri');
  
  mesh.TriC           = ncread(filename,'TriC');
  mesh.Tricc          = ncread(filename,'Tricc');
  
  mesh.VAc            = ncread(filename,'VAc');
  mesh.nAc            = size( mesh.VAc,1);
  mesh.Aci            = ncread(filename,'Aci');
  mesh.iAci           = ncread(filename,'iAci');
  
  mesh.A              = ncread(filename,'A');
  mesh.R              = ncread(filename,'R');

end