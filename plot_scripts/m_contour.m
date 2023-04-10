addpath('~/jwplot/common/')

load1 = load('../matfile.mat');

matf = load1.matfe;
lsb = load1.lsb;
lsh = load1.lsh;

cellM = cell(numel(lsh),1);

matf = matf';
for bi = 1:numel(lsb)
    fe = matf(:,bi);
    [lsh1, z1] = calM_fe2M(lsh,fe,1);
    cellM{bi} = reshape(z1,1,[]);
end

zm = cell2mat(cellM);
zm = zm';
x = lsh1;
y = lsb.^-1;

xq = x;
yq = y;

[xm, ym] = ndgrid(x,y);
[xmq, ymq] = ndgrid(xq,yq);
zmq = interpn(xm,ym,zm,xmq,ymq,'spline');

f = figure;
contourf(xmq,ymq,zmq,32,'LineStyle','none')

ax = gca;

clb = colorbar;
title(clb,'$M$','Interpreter','latex')

box on
grid off
view(2)
shading interp
colormap turbo

set(gca,'YScale','log')

f.Position(3:4) = 420.*[1,1];
pbaspect([1 1 1])

xlabel('$h/J$','Interpreter','latex')
ylabel('$T/J$','Interpreter','latex')

ax.FontSize = 16;