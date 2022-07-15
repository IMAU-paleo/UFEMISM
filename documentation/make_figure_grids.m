clc
clear all
close all

wa = 400;
ha = 400;

fs = 18;

H.Fig = figure('color','w','position',[500,500,wa+75,ha]);
H.Ax  = axes('units','pixels','position',[0,0,wa,ha],'xtick',[],'ytick',[],'xlim',[-1.3,1.3],'ylim',[-1.3,1.3]);

% cell i,j
x = [-0.5,0.5,0.5,-0.5,-0.5];
y = [-0.5,-0.5,0.5,0.5,-0.5];
line('xdata',x,'ydata',y,'linestyle','--');
V = [-0.5,-0.5;0.5,-0.5;0.5,0.5;-0.5,0.5];
f = [1,2,3,4];
patch('faces',f,'vertices',V,'facecolor','k','facealpha',0.1,'edgecolor','none');
for xi = -1:1
  y = -1:1;
  x = [0,0,0]+xi;
  line('xdata',x,'ydata',y,'linestyle','-');
end
for yj = -1:1
  x = -1:1;
  y = [0,0,0]+yj;
  line('xdata',x,'ydata',y,'linestyle','-','color','k','marker','o','markerfacecolor','k','markeredgecolor','k','markersize',12);
end


istring = {'i-1','i','i+1'};
jstring = {'j-1','j','j+1'};

% A
for i=-1:1
  for j=-1:1
    text(i+0.04,j+0.1,['A(' istring{i+2} ',' jstring{j+2} ')'],'interpreter','latex','fontsize',fs);
  end
end

% B
for i=-0.5:0.5
  for j = -0.5:0.5
    line('xdata',i,'ydata',j,'linestyle','none','marker','o','markerfacecolor','b','markeredgecolor','k','markersize',12);
  end
end
text(0.55,0.6,'B(i,j)','interpreter','latex','fontsize',fs);

% Cx
for i=-0.5:0.5
  for j = -1:1
    line('xdata',i,'ydata',j,'linestyle','none','marker','o','markerfacecolor','r','markeredgecolor','k','markersize',12);
  end
end
text(0.55,0.1,'Cx(i,j)','interpreter','latex','fontsize',fs);

% Cy
for i=-1:1
  for j = -0.5:0.5
    line('xdata',i,'ydata',j,'linestyle','none','marker','o','markerfacecolor',[0,0.6,0],'markeredgecolor','k','markersize',12);
  end
end
text(0.05,0.6,'Cy(i,j)','interpreter','latex','fontsize',fs);