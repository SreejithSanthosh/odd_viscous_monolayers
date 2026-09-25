function c = coeff(location,~)
% This function defines the coefficients used to solve the 
% constitutive relation with spatiotemporally varying nematic order 

Nz = numel(location.x);
xArr = location.x; yArr = location.y;
[s,phi] = orient_order(xArr,yArr); 

c = zeros([16,Nz]); [material_moduli,~] = get_parameters;

bulk_visc = material_moduli(1); shear_visc = material_moduli(2);lambda_0 = material_moduli(3);
lambda_s = material_moduli(4); del_nu = material_moduli(5);

% Adding the order parameter dependence and spatial array
bulk_visc = bulk_visc+0*s; shear_visc = shear_visc+0*s; lambda_0 = lambda_0*s;
lambda_s = lambda_s*s; del_nu = del_nu*s;

%C11
c(1,:) = (bulk_visc+shear_visc+del_nu.*cos(4*phi))/2+lambda_s.*cos(2*phi);
c(2,:) = (0.5*lambda_s+0.5*lambda_0+del_nu.*cos(2*phi)).*sin(2*phi);
c(3,:) = (0.5*lambda_s-0.5*lambda_0+del_nu.*cos(2*phi)).*sin(2*phi);
c(4,:) = 0.5*shear_visc-0.5*del_nu.*cos(4*phi);

%C21
c(5,:) = (0.5*lambda_s+0.5*lambda_0+del_nu.*cos(2*phi)).*sin(2*phi);
c(6,:) = 0.5*bulk_visc-0.5*shear_visc-lambda_0.*cos(2*phi)-0.5*del_nu.*cos(4*phi);
c(7,:) = 0.5*shear_visc-0.5*del_nu.*cos(4*phi);
c(8,:) = (0.5*lambda_s-0.5*lambda_0-del_nu.*cos(2*phi)).*sin(2*phi);

%C12
c(9,:) = (0.5*lambda_s-0.5*lambda_0+del_nu.*cos(2*phi)).*sin(2*phi);
c(10,:) = 0.5*shear_visc-0.5*del_nu.*cos(4*phi);
c(11,:) = 0.5*bulk_visc-0.5*shear_visc+lambda_0.*cos(2*phi)-0.5*del_nu.*cos(4*phi);
c(12,:) = (0.5*lambda_s+0.5*lambda_0-del_nu.*cos(2*phi)).*sin(2*phi);

%C22
c(13,:) = 0.5*shear_visc-0.5*del_nu.*cos(4*phi);
c(14,:) = (0.5*lambda_s-0.5*lambda_0-del_nu.*cos(2*phi)).*sin(2*phi);
c(15,:) = (0.5*lambda_s+0.5*lambda_0-del_nu.*cos(2*phi)).*sin(2*phi);
c(16,:) = (bulk_visc+shear_visc+del_nu.*cos(4*phi))/2-lambda_s.*cos(2*phi);

end 