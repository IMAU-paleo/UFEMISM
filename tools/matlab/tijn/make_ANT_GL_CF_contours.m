clc
close all
clear all

filename = '/Users/berends/Documents/Datasets/BedMachine_Antarctica/Bedmachine_v1_Antarctica_2km.nc';

filename_dst_GL = 'ANT_GL_hires.mat';
filename_dst_CF = 'ANT_CF_hires.mat';

x  = ncread( filename,'x');
y  = ncread( filename,'y');
Hi = ncread( filename,'Hi');
Hb = ncread( filename,'Hb');

% Calculate thickness above flotation
seawater_density = 1014;
ice_density      = 910;
TAF = Hi - max( 0, -Hb * seawater_density / ice_density);

M_GL = process_contour( contour( x,y,TAF',[0,0]));
GL.x = M_GL( 1,:);
GL.y = M_GL( 2,:);
M_CF = process_contour( contour( x,y,Hi' ,[0.1,0.1]));
CF.x = M_CF( 1,:);
CF.y = M_CF( 2,:);
close all

if exist( filename_dst_GL,'file')
  delete( filename_dst_GL)
end
save( filename_dst_GL,'GL');

if exist( filename_dst_CF,'file')
  delete( filename_dst_CF)
end
save( filename_dst_CF,'CF');

function M2 = process_contour( M)

M2 = [];
ii = 1;
while ii < size( M,2)
  n = M(2,ii);
  M2 = [M2, M(:,ii+1:ii+n), [NaN;NaN]];
  ii = ii+n+1;
end

end