function [vec3] =blochPlot3(vec,col,hon,fig,yo) 

% init = 0;
% 
% vec=[1 0.72*1i];
% col=[0.4 0.3 0.9];
% 
% hon = 0; 
% fig=3;

if ~exist('yo')
    yo=1;
end

myaapub = 1;

nnn=0;
% yo=2;

vec = vec / sqrt(vec*vec');
rho = vec'*vec;

sigx = [0 1;1 0];
sigy = [0 -1i;1i 0];
sigz = [ 1 0;0 -1];

x = trace(rho*sigx);
y = trace(rho*sigy);
z = trace(rho*sigz);

vec3 = [x y z];


a=0;
for p=linspace(0,2*pi,500)
    a=a+1;
    x4(a) = cos(p);
    y4(a) = sin(p);
    z4(a) = 0;
    x5(a) = 0;
    y5(a) = sin(p);
    z5(a) = cos(p);
    x7(a) = sin(p);
    y7(a) = 0;
    z7(a) = cos(p);
end

lw=2;
h1=figure(fig); 
set(0,'defaultaxesposition',[0 0 1 1])
figpos = [25 25 600 600];
set(gcf, 'Position', figpos);


if hon==0
hold('off');
sphere3(30); 
colormap winter; 
hold on
plot3(0,0,0,'.','MarkerSize',15*yo,'Color',[0 0 0])

plot3(x4,y4,z4,'-','LineWidth',2,'Color',[0 0 0])
plot3(x5,y5,z5,'-','LineWidth',2,'Color',[0 0 0])
plot3(x7,y7,z7,'-','LineWidth',2,'Color',[0 0 0])
plot3([-1 1],[0 0] ,[0 0] ,'-','LineWidth',2,'Color',[0 0 0])
plot3([0 0],[-1 1] ,[0 0] ,'-','LineWidth',2,'Color',[0 0 0])
plot3([0 0],[0 0] ,[-1 1] ,'-','LineWidth',2,'Color',[0 0 0])
grid off
axis image
axis off
axis vis3d
set(gcf,'color','w');
alpha(0.5); 
light('Position',[1 -1 1]);
else
    hold on;
end

plot3([0 x],[0 y],[0 z],'-','Color',col,'LineWidth',2)
 plot3(x,y,z,'.','MarkerSize',15*yo,'Color',col)


view([45 35])

% pause(4)
% if myaapub == 1
%     
% %      axis image
% %     set(gcf, 'Position', figpos);
%  myaa('publish')
% 
% %     print(gcf,'-dpng','illmub');
%     
% %     figure(3)
% %     myaa('publish');
% %     figure(1)
% %     f=getframe(gcf);
% %     imwrite(f.cdata,[ num2str(nnn) '.png']);
% else
%     print(gcf,'-dpng',['/Users/eliotbolduc/Desktop/Gmovie/bloch']);
end


