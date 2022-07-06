function d_trans = ReadTransectData( filename, varname, ti)

vi_transect    = ncread( filename,'vi_transect');
nV_transect    = size( vi_transect,1);
w_transect     = ncread( filename,'w_transect');

d_trans = zeros( nV_transect,1);
d = ncread( filename, varname, [1,ti],[Inf,1]);

for vii = 1: nV_transect
  vi = vi_transect( vii,1);
  vj = vi_transect( vii,2);
  wi = w_transect(  vii,1);
  wj = w_transect(  vii,2);
  
  d_trans( vii) = d(vi)*wi + d(vj)*wj;
end

end