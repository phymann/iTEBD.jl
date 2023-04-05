addpath('~/jwplot/common')

load1 = load('../matfile.mat');

matf = load1.matfe;
lsb = load1.lsb;
lsh = load1.lsh;

cellcv = cell(numel(lsh),1);
lsmax = lsh;

matf = matf';
for hi = 1:numel(lsh)
    fe = matf(hi,:);
    
    [lst, z1] = calCv_fe2S2Cv(lsb,fe,1);

    lst_sp = linspace(lst(1),lst(end),3*numel(lst));
    z1_sp = interp1(lst,z1,lst_sp,'spline');
    [~,idx] = max(z1_sp);
    lsmax(hi) = lst_sp(idx);
    
    cellcv{hi} = reshape(z1_sp,1,[]);
end

zm = cell2mat(cellcv);
x = lsh;
y = lst_sp;

xq = x;
yq = y;

[xm, ym] = ndgrid(x,y);
[xmq, ymq] = ndgrid(xq,yq);
zmq = interpn(xm,ym,zm,xmq,ymq,'spline');

f = figure;
contourf(xmq,ymq,zmq,256,'LineStyle','none')

hold on

plot(lsh,lsmax,'.','Color','k')

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

% find local maximal for fixed beta
x_sp = linspace(x(1),x(end),2*numel(x));
lshmax = y;
for i = 1:size(zmq,2)
    z_sp = interp1(x,zm(:,i),x_sp,'spline');
    [~,idx] = max(z_sp);
    lshmax(i) = x_sp(idx);
end
% hold on
% plot(lshmax,y,'*','Color','b')

clim([0, 1])
