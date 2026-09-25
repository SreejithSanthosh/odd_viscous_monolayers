function bc = bc_sfbc(location,~)
% The stress free boundary condition 
[~,~,params_sfbc] = get_parameters;

% SFBC from isotropic prestress
m_iso = params_sfbc(1); m_aniso = params_sfbc(2);
bc_iso = -m_iso*[location.nx;location.ny];

% Load nematic field at boundary 
xArr = location.x; yArr = location.y;
[s,phi] = orient_order(xArr,yArr);

% SFBC from anisotropic prestress
q1 = s.*cos(2*phi)*0.5; q2 = s.*sin(2*phi)*0.5;
bc_aniso = -m_aniso*[q1.*location.nx+q2.*location.ny;...
    q2.*location.nx-q1.*location.ny];

bc = bc_iso+bc_aniso;

end 


