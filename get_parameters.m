function [material_moduli,kappa,m_iso] = get_parameters

% Material constant
bulk_visc = 1; shear_visc = 10^(0); 
lambda_0 = 0.0; lambda_s = 0.2; 

del_nu = 0;
material_moduli = [bulk_visc,shear_visc,lambda_0,lambda_s,del_nu];

% Boundary condition 
kappa = 10^(-4); 

% Pre-stress coefficient 
m_iso = 1; 
end 
