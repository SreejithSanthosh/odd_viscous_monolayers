function f = activeforcing(location,~)

Nz = numel(location.x); xArr = location.x; yArr = location.y;

dxy = 10^(-4);
[~,~,params_prestress] = get_parameters;
m_aniso = params_prestress(2);


[s,phi] = orient_order(xArr,yArr);
[s_xp,phi_xp] = orient_order(xArr+dxy,yArr);
[s_yp,phi_yp] = orient_order(xArr,yArr+dxy);

alpha1_x = 0*s; alpha1_y = 0*s;
alpha2_x = 0.5*m_aniso*(s_xp.*cos(2*phi_xp)-s.*cos(2*phi))./dxy; 
alpha2_y = 0.5*m_aniso*(s_yp.*cos(2*phi_yp)-s.*cos(2*phi))./dxy;

alpha3_x = 0.5*m_aniso*(s_xp.*sin(2*phi_xp)-s.*sin(2*phi))./dxy; 
alpha3_y = 0.5*m_aniso*(s_yp.*sin(2*phi_yp)-s.*sin(2*phi))./dxy;

f = zeros([2,Nz]); 
f(1,:) = alpha3_y+alpha1_x+alpha2_x; 
f(2,:) = alpha3_x+alpha1_y-alpha2_y;

end 