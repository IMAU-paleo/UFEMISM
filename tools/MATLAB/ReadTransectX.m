function x = ReadTransectX( filename)

vi_transect    = ncread( filename,'vi_transect');
nV_transect    = size( vi_transect,1);
w_transect     = ncread( filename,'w_transect');

x = zeros( nV_transect,1);
V = ncread( filename, 'V');

for vii = 1: nV_transect
  vi = vi_transect( vii,1);
  vj = vi_transect( vii,2);
  wi = w_transect(  vii,1);
  wj = w_transect(  vii,2);
  
  x( vii) = V(vi,1)*wi + V(vj,1)*wj;
end

end