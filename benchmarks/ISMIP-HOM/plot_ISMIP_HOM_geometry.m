clc
clear all
close all

L = 160*1e3;

x = linspace( -L/2, L/2, 30);
y = x;

[yg,xg] = meshgrid(y,x);

Hs = 2000.0 - xg * tan( 0.5 * pi / 180.0);
Hb = Hs - 1000.0 + 500.0 * sin( xg * 2.0 * pi / L) .* sin( yg * 2.0 * pi / L);
Hi = Hs - Hb;

surf(Hb,'facecolor',[0.8,0.4,0.1]);
hold on
surf(Hs,'facecolor','b','facealpha',0.5);
set(gca,'cameraposition',[0.3356    0.3237    4.0721]*1e3,'xtick',[],'ytick',[],'ztick',[]);
set(gcf,'position',[744   370   758   680],'color','w');
