% FEM for odd elastic solid - no time dependence

clear; clc; close all;

% We are making an 
model = createpde(2); % 2 is the number of equatiosn in the PDE 
geometryFromEdges(model,@circleg);

% simulation of the material 
specifyCoefficients(model,"m",0,"d",0,"c",@coeff,"a",0,"f",@activeforcing);
applyBoundaryCondition(model,"neumann","edge",1:model.Geometry.NumEdges,"q",@bc1,"g",@bc2);

clc; hmax = 0.01; generateMesh(model,"Hmax",hmax);
results = solvepde(model);

clc; N = 100; xq = linspace(-1,1,N); yq = xq;
[X,Y] = ndgrid(xq,yq); xArr = X(:); yArr = Y(:);

uinterp = interpolateSolution(results,xArr,yArr,1);
vinterp = interpolateSolution(results,xArr,yArr,2);

uinterp = reshape(uinterp,size(X)); 
vinterp = reshape(vinterp,size(X));

[ux,uy] = evaluateGradient(results,xArr,yArr,1);
[vx,vy] = evaluateGradient(results,xArr,yArr,2);

S11 = 0.5*(ux-vy); S21 = 0.5*(uy+vx);
%Compute the angles
EdMag = 0*S11; EdPhi = EdMag;
for i = 1:N^2
    E11 = S11(i); E12 = S21(i);
    Ed = [E11,E12; E12,-E11];
    [V, D] = eig(Ed); [~,ind] = sort(diag(D));
    Vs = V(:,ind); Phi = atan2(Vs(2,end),Vs(1,end));
    EdMag(i) = max(diag(D)); EdPhi(i) = Phi;
end


%% Visualize the results
close all; clc
% Draw the circles
phiCirc = 0:0.001:2*pi; phiCirc = [phiCirc,2*pi]; rCirc = 0.99;
xCirc = rCirc*cos(phiCirc); yCirc = rCirc*sin(phiCirc);

% Coarsen the quiver plots 
[S,phi] = orient_order(X,Y);
S(sqrt(X.^2+Y.^2)>1) = nan; phi(sqrt(X.^2+Y.^2)>0.92) = nan;

coarse = 7;
XCr = X(1:coarse:end,1:coarse:end);YCr = Y(1:coarse:end,1:coarse:end);
nxCr = cos(phi(1:coarse:end,1:coarse:end)); nyCr = sin(phi(1:coarse:end,1:coarse:end));
v1Cr = uinterp(1:coarse:end,1:coarse:end); v2Cr = vinterp(1:coarse:end,1:coarse:end);

axlim = 1.1; % This is the coordinate limits;
fntSz = 24; % Font Size
linSz = 2;

figure('Position',[15 625 2500 550],'color','w')

subplot(1,4,1)
contourf(X,Y,S,'EdgeColor','none'); axis equal off
hold on; plot(xCirc,yCirc,'k','LineWidth',linSz+1)
quiver(XCr,YCr,nxCr,nyCr,'m','Linewidth',linSz,'ShowArrowHead','off','AutoScaleFactor',0.7); hold off
title('$ \mathbf{Q}(\mathbf{x}) $','Interpreter','latex','FontSize',fntSz,'Color','k'); 
c = colorbar; c.FontSize = fntSz; c.LineWidth = linSz; c.Color = 'k';
xlim([-axlim axlim]); ylim([-axlim axlim]);

% Changing colorbar properties
c.Position(1) = c.Position(1) +0.015;
c.Position(2) = c.Position(2) - 0.01;
c.Position(4) = c.Position(4) * 0.75;
c.Position(3) = c.Position(3) * 1.25;
clim([0,1])



subplot(1,4,2)
vMag = sqrt(uinterp.^2+vinterp.^2);
contourf(X,Y,vMag,'EdgeColor','none'); axis equal off
hold on; plot(xCirc,yCirc,'k','LineWidth',linSz+1)
quiver(XCr,YCr,v1Cr,v2Cr,'k','Linewidth',linSz,'AutoScaleFactor',0.7); hold off
title('$ \mathbf{v}(x) $','Interpreter','latex','FontSize',fntSz,'Color','k'); 
c = colorbar; c.FontSize = fntSz; c.LineWidth = linSz; c.Color = 'k';
xlim([-axlim axlim]); ylim([-axlim axlim]);

% Changing colorbar properties
c.Position(1) = c.Position(1) +0.015;
c.Position(2) = c.Position(2) - 0.01;
c.Position(4) = c.Position(4) * 0.75;
c.Position(3) = c.Position(3) * 1.25;


isoStrain = reshape(ux+vy,[N,N]); 
ax = subplot(1,4,3);
ax.FontSize = fntSz;
contourf(X,Y,isoStrain,'EdgeColor','none'); axis equal off
hold on; plot(xCirc,yCirc,'k','LineWidth',linSz+1)
title('Isotropic strain','Interpreter','latex','FontSize',fntSz,'Color','k'); c = colorbar;
c.FontSize = fntSz; c.LineWidth = linSz; c.Color = 'k';
xlim([-axlim axlim]); ylim([-axlim axlim]);


% Changing colorbar properties
c.Position(1) = c.Position(1) +0.015;
c.Position(2) = c.Position(2) - 0.01;
c.Position(4) = c.Position(4) * 0.75;
c.Position(3) = c.Position(3) * 1.25;

% Compute stress 
beta1 = (ux+vy)*0.5; beta2 = (ux-vy)*0.5; beta3 = (uy+vx)*0.5;
beta1 = reshape(beta1,[N,N]); beta2 = reshape(beta2,[N,N]); beta3 = reshape(beta3,[N,N]);
[a1,a2,a3] = compute_stress(beta1,beta2,beta3,X,Y);
[paramsMaterial,paramsBoundary,params_prestress] = get_parameters;

ax = subplot(1,4,4);
ax.FontSize = fntSz;
isoStress = reshape(a1,[N,N])+params_prestress(1);
contourf(X,Y,isoStress,'EdgeColor','none'); axis equal off
hold on; plot(xCirc,yCirc,'k','LineWidth',linSz+1)
title('Isotropic stress','Interpreter','latex','FontSize',fntSz,'Color','k')
c = colorbar; c.Color = 'k';
c.FontSize = fntSz; c.LineWidth = linSz;
% 
% cMax = max(abs(isoStress(:)));
% clim([-cMax cMax])
xlim([-axlim axlim]); ylim([-axlim axlim]);

% Changing colorbar properties
c.Position(1) = c.Position(1) +0.015;
c.Position(2) = c.Position(2) - 0.01;
c.Position(4) = c.Position(4) * 0.75;
c.Position(3) = c.Position(3) * 1.25;

sgtitle(sprintf('$ \\lambda^s = %.1f $ $ \\lambda^o = %.1f $',...
     paramsMaterial(4),paramsMaterial(3)),'Interpreter','latex','FontSize',fntSz,'Color','k');

exportgraphics(gcf,'plot_fig4paper.png','Resolution',600)