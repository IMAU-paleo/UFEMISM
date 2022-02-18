function A = CSR_to_sparse( m, n, A_ptr, A_index, A_val)

nnz = A_ptr(end)-1;

Ai = zeros( nnz, 1);
Aj = zeros( nnz, 1);
Av = zeros( nnz, 1);

k = 0;
for i = 1: m
  for ii = A_ptr( i): A_ptr( i+1)-1
    j = A_index( ii);
    v = A_val(   ii);
    k = k+1;
    Ai( k) = i;
    Aj( k) = j;
    Av( k) = v;
  end
end

A = sparse(Ai,Aj,Av,m,n);

end