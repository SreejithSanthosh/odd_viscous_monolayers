function bcCond = bc_drag(location,~)
% The dirichlet componenet of boundary condition to generate 
% generate a finite velocity field (constrain the null space of possible velocity fields)
nx = location.nx; ny = location.ny;
[~,kappa,~] = get_parameters;
bcCond = zeros(2,2);

bcCond(1,1) = kappa*nx.^2+kappa*(1-nx.^2);
bcCond(1,2) = kappa*nx.*ny-kappa.*nx.*ny;
bcCond(2,1) = kappa*nx.*ny-kappa.*nx.*ny;
bcCond(2,2) = kappa*ny.^2+kappa*(1-ny.^2);

end 

