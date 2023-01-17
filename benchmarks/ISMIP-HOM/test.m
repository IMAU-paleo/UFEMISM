clc
clear all
close all


H.Fig = figure;
H.Ax  = axes;

x = linspace(0,2*pi,100);
y = linspace(0,2*pi,100);
[xg,yg] = meshgrid(x,y);
c = sin(xg).*sin(yg);

H.im = imagesc('xdata',x,'ydata',y,'cdata',c);
set(gca,'clim',[-1,1]);

for t = linspace(0,10*pi,1000)
  c = sin(t)*sin(xg).*sin(yg);
  set(H.im,'cdata',c);
  drawnow('update');
  pause(0.1)
end