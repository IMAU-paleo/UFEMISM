function average_f = vertical_average( zeta, f)
% Calculate the vertical average of any given function f defined at the vertical zeta grid.

average_f = 0;

nz = length( zeta);

for k = 1: nz-1
  average_f = average_f + 0.5 * (f(k+1) + f(k)) * (zeta(k+1) - zeta(k));
end

end