function [material_moduli,kappa,params_sfbc] = get_parameters

% Material constant
bulk_visc = 1; shear_visc = 10^(0); 
lambda_0 = 0; lambda_s = 0.0; 

del_nu = 0;
material_moduli = [bulk_visc,shear_visc,lambda_0,lambda_s,del_nu];

% Boundary condition 
kappa = 10^(-4); 

% Pre-stress coefficient 
m_iso = 0; m_aniso = 1;
params_sfbc = [m_iso, m_aniso];
end 
