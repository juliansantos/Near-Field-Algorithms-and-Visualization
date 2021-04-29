% This function implement a initial guess for the phase of the field
% desired. 
function initialguess(x_mesh, y_mesh, M, layers)
    ap_dim =[7.2,10.13]; % Size of the aperture of the aut
    
    % Initial guess according to Razavi_2007
    m1 = -ap_dim(1)/2<=x_mesh & x_mesh<=ap_dim(1)/2;
    m2 = -ap_dim(2)/2<=y_mesh & y_mesh<=ap_dim(2)/2;
    Maut = m1 & m2;
    
    % Propagation of Maut to layer 1
    
    
end