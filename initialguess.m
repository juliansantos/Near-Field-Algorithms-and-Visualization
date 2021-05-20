% This function implement a initial guess for the phase of the field
% desired. 
function Maut = initialguess(x_mesh, y_mesh, field, layer_1)
    ap_dim =[7.2,10.13]; % Size of the aperture of the aut
        
    % Initial guess: (Basic approach) 1<0 iff x,y E aperture otherwise 0
    % Method = Basic Approach
    
    Method = 2;
        % Method = 1 -> Basic approach
        % Method = 2 -> Hybrid approach (requires global optimization technique)
        % Method = 3 -> Binary magnitude image
    
    m1 = -ap_dim(1)/2<=x_mesh & x_mesh<=ap_dim(1)/2;
    m2 = -ap_dim(2)/2<=y_mesh & y_mesh<=ap_dim(2)/2;
    Maut = m1 & m2;
       
    if Method == 2 
        lambda = 0.005; % TODO; put as argument
        k0 = 2*pi/lambda; % wave number 
        theta = pi/3;% Theta estimated for the initial guess
        phi = pi/4;% Phi estimated for the initial guess
        
        Maut = Maut .* exp(-1i*k0*(x_mesh*sin(theta)*cos(phi) + y_mesh*sin(theta)*sin(phi)));
    elseif Method == 3
       bmi = im2bw(abs(f_mesh(:,:,1))/max(max(abs(f_mesh(:,:,1)))),1/sqrt(2));
       x_center = mean(x_mesh(bmi==1)); % Orientation center
       y_center = mean(y_mesh(bmi==1)); % Orientation center
       [TH,PHI,~] = cart2sph(x_center,y_center,layer_1) ;
       I = exp(-1i*k0*(x_mesh*sin(TH)*cos(PHI) + y_mesh*sin(TH)*sin(PHI)));
       
       % convert current to electric field. 
       
    end
    % Initial guess: (evolutive algorithm approach)
    
end