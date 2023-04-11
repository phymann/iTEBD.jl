logornot = 2;

addpath('~/jwplot/common/')

% load1 = load('../matfile.mat');
load1 = load('../D20230411T102558.mat');

matf = load1.matfe;
lsb = load1.lsb;
lsh = load1.lsh;

cellM = cell(numel(lsh),1);
lsmax = lsh;

matf = matf';
for bi = 1:numel(lsb)
    fe = matf(:,bi);
    [lsh1, z1] = calchim_fe2chim(lsh,fe,1,1,'spline');

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
switch logornot
    case 1
        contourf(xmq,ymq,log(zmq),256,'LineStyle','none')
    case 2
        contourf(xmq,ymq,(zmq),256,'LineStyle','none')
end

ax = gca;

clb = colorbar;
title(clb,'$\chi_m$','Interpreter','latex')

box on
grid off
view(2)
shading interp
colormap hot

% set(gca,'YScale','log')
% caxis([0,10])

f.Position(3:4) = 420.*[1,1];
pbaspect([1 1 1])

xlabel('$h/J$','Interpreter','latex')
ylabel('$T/J$','Interpreter','latex')

ax.FontSize = 16;

% find local maximal for fixed h
y_sp = linspace(y(1),y(end),2*numel(y));
lstmax = x;
for i = 1:size(zmq,1)
    z_sp = interp1(y,zm(i,:),y_sp,'spline');
    [~,idx] = max(z_sp);
    lstmax(i) = y_sp(idx);
end
hold on
plot(x,lstmax,'.','Color','w')

% find local maximal for fixed beta
x_sp = linspace(x(1),x(end),2*numel(x));
lshmax = y;
for i = 1:size(zmq,2)
    z_sp = interp1(x,zm(:,i),x_sp,'spline');
    [~,idx] = max(z_sp);
    lshmax(i) = x_sp(idx);
end
hold on
plot(lshmax,y,'*','Color','g')

text(-0.48, 2.5, 'max for fixed \beta', 'Color', 'g', 'FontSize', 15)
text(-0.48, 2.8, 'max for fixed h', 'Color', 'w', 'FontSize', 15)

clim([0,100])
