% Author : Sreejith Santhosh 
% Created: 9/21/26
% The parameters are defined in get_parameters.m: Use that file to change
% the material properties 
% To change the nematic field being probed, change it using the file
% orient_order.m

clear; clc; close all;
fntSz = 20; linSz = 1;

%% Using PDE toolbox 
% Define Eqns and geometry 
model = createpde(2); 
geometryFromEdges(model,@circleg);

% Specify the coefficients 
specifyCoefficients(model,"m",0,"d",0,"c",@coeff,"a",0,"f",[0,0]');
applyBoundaryCondition(model,"neumann","edge",1:model.Geometry.NumEdges,"q",@bc_drag,"g",@bc_sfbc);

% Generate mesh and solve the equation 
hmax = 0.01; generateMesh(model,"Hmax",hmax);
results = solvepde(model);

%% Visualize the results
clc; N = 100; xq = linspace(-1,1,N); yq = xq;
[X,Y] = ndgrid(xq,yq); xArr = X(:); yArr = Y(:);

uinterp = interpolateSolution(results,xArr,yArr,1);
vinterp = interpolateSolution(results,xArr,yArr,2);

uinterp = reshape(uinterp,size(X)); 
vinterp = reshape(vinterp,size(X));

[ux,uy] = evaluateGradient(results,xArr,yArr,1);
[vx,vy] = evaluateGradient(results,xArr,yArr,2);

% Principal magnitude and axis of the deviatoric strain-rate tensor.
S11 = 0.5 * (ux - vy);
S21 = 0.5 * (uy + vx);
EdMag = hypot(S11, S21);
EdPhi = 0.5 * atan2(S21, S11);
EdPhi(EdMag == 0) = NaN; % No unique principal axis at zero magnitude

[fig, ax] = plot_fields(X, Y, uinterp, vinterp, ux, uy, vx, vy, ...
    EdMag, EdPhi, fntSz, linSz);
