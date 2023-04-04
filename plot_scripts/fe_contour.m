addpath('~/jwplot/common')

load1 = load('../matfile.mat');

matf = load1.matfe;
lsb = load1.lsb;
lsh = load1.lsh;

cellcv = cell(numel(lsh),1);
lsmax = lsh;

matf = matf';

zm = matf;
x = lsh;
y = lsb.^-1;

xq = x;
yq = y;

[xm, ym] = ndgrid(x,y);
[xmq, ymq] = ndgrid(xq,yq);
zmq = interpn(xm,ym,zm,xmq,ymq,'spline');

f = figure;
contourf(xmq,ymq,zmq,256,'LineStyle','none')


ax = gca;

clb = colorbar;
title(clb,'$C_m$','Interpreter','latex')

box on
grid off
view(2)
shading interp
colormap hot

f.Position(3:4) = 420.*[1,1];
pbaspect([1 1 1])

xlabel('$h/J$','Interpreter','latex')
ylabel('$T/J$','Interpreter','latex')

ax.FontSize = 16;
