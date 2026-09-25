function bc_iso = bc_sfbc(location,~)
% The stress free boundary condition 
[~,~,m_iso] = get_parameters;

% SFBC from isotropic prestress
bc_iso = -m_iso*[location.nx;location.ny];

end 


