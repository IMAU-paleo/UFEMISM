function A = read_CSR_from_NetCDF( filename)

% Get matrix size
f = ncinfo( filename);
m   = [];
n   = [];
nnz = [];
for di = 1: length( f.Dimensions)
  if     strcmpi(f.Dimensions(di).Name,'m')
    m   = f.Dimensions(di).Length;
  elseif strcmpi(f.Dimensions(di).Name,'n')
    n   = f.Dimensions(di).Length;
  elseif strcmpi(f.Dimensions(di).Name,'nnz')
    nnz = f.Dimensions(di).Length;
  end
end

% Safety
if isempty(m) || isempty(n) || isempty(nnz)
  error('Couldnt find matrix size in file!')
end

% Read matrix data
ptr   = ncread( filename, 'ptr'  );
index = ncread( filename, 'index');
val   = ncread( filename, 'val'  );

% Convert to Matlab sparse matrix
A = CSR_to_sparse( m, n, ptr, index, val);

end