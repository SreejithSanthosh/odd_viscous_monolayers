function [alpha1,alpha2,alpha3] = compute_stress(beta1,beta2,beta3,X,Y)
[s,phi] = orient_order(X,Y); 

[paramsMaterial,~] = get_parameters;

B = paramsMaterial(1); mu = paramsMaterial(2);lambA = paramsMaterial(3);
lamb0 = paramsMaterial(4); del_mu = paramsMaterial(5);

%Computing the akoha components
alpha1 = B*beta1+(lamb0-lambA).*s.*cos(2*phi).*beta2+(lamb0-lambA).*s.*sin(2*phi).*beta3;
alpha2 = (lamb0+lambA).*s.*cos(2*phi).*beta1+(mu+s.*del_mu.*cos(4*phi)).*beta2+s.*del_mu.*sin(4*phi).*beta3;
alpha3 = (lamb0+lambA).*s.*sin(2*phi).*beta1+s.*del_mu.*sin(4*phi).*beta2+(mu-s.*del_mu.*cos(4*phi)).*beta3;


end 