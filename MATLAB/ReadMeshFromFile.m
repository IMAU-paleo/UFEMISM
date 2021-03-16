function mesh = ReadMeshFromFile(filename)

  mesh.V              = ncread(filename,'V');
  mesh.nV             = size(mesh.V,1);
  mesh.Tri            = ncread(filename,'Tri');
  mesh.nTri           = size(mesh.Tri,1);
  
  mesh.xmin           = mesh.V(1,1);
  mesh.xmax           = mesh.V(2,1);
  mesh.ymin           = mesh.V(2,2);
  mesh.ymax           = mesh.V(3,2);
  
  mesh.nC             = ncread(filename,'nC');
  mesh.C              = ncread(filename,'C');
  mesh.R              = ncread(filename,'R');
  
  mesh.edge_index     = ncread(filename,'edge_index');
  mesh.Tri_edge_index = ncread(filename,'Tri_edge_index');
  
  mesh.niTri          = ncread(filename,'niTri');
  mesh.iTri           = ncread(filename,'iTri');
  
  mesh.TriC           = ncread(filename,'TriC');
  
  mesh.vi_transect    = ncread(filename,'vi_transect');
  mesh.nV_transect    = size(mesh.vi_transect,1);
  mesh.w_transect     = ncread(filename,'w_transect');

end