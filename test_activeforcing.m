% This code is to test the activeforcing function 

clc; hmax = 0.05; generateMesh(model,"Hmax",hmax);
results = solvepde(model);

clc; N = 200; xq = linspace(-1,1,N); yq = xq;
[X,Y] = ndgrid(xq,yq); xArr = X(:); yArr = Y(:);

location.x = xArr; location.y = yArr;
f = activeforcing(location,location);

fx =reshape(f(1,:),[N,N]);
fy =reshape(f(2,:),[N,N]);

fx(sqrt(X.^2+Y.^2)>1) = nan;
fy(sqrt(X.^2+Y.^2)>1) = nan;

cr = 10;
contourf(X,Y,sqrt(fx.^2+fy.^2),'edgecolor','none'); hold on 
quiver(X(1:cr:end,1:cr:end),Y(1:cr:end,1:cr:end),...
    fx(1:cr:end,1:cr:end),fy(1:cr:end,1:cr:end),'r'); hold off; axis equal off
c  = colorbar; c.FontSize = 20; set(gcf,'color','w')