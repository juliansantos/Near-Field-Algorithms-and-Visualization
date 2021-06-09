% This function implement a initial guess for the phase of the field
% desired. 
function [Faut, layer_aut] = initialguess(X, Y, Z, Ex, Method, layers,plot)
    ap_dim =[7.2,10.13]; % Size of the aperture of the aut
    zref = 41.2; % Relative position of the center of the aperture in z   
    layer_aut = 0;
    lambda = 0.005; 
    k0 = 2*pi/lambda;
    % Initial guess: (Basic approach) 1<0 iff x,y E aperture otherwise 0
    % Method = Basic Approach
    
    %Method = 1;
        % Method = 0 -> No Initial Guess
        % Method = 1 -> Basic approach
        % Method = 2 -> Nearest Layer to AUT in Data (Magnitude)
        % Method = 3 -> Hybrid approach (requires global optimization technique)
        % Method = 4 -> Binary magnitude image
    [x_mesh, y_mesh, ~]= getFieldLayer(X, Y, Z, Ex, 1);

    if Method == 1
        m1 = -ap_dim(1)/2<=x_mesh & x_mesh<=ap_dim(1)/2;
        m2 = -ap_dim(2)/2<=y_mesh & y_mesh<=ap_dim(2)/2;
        Faut = m1 & m2; 
    elseif Method == 2 
        [~,idx]=min(abs(0 + zref - Z)); layer_aut= Z(idx);
        [~, ~, Faut]= getFieldLayer(X, Y, Z, Ex, layer_aut-zref);
        Faut = abs(Faut);
        layer_aut = layer_aut - zref;
    elseif Method == 3
        lambda = 0.005; % TODO; put as argument
        k0 = 2*pi/lambda; % wave number 
        theta = pi/3;% Theta estimated for the initial guess
        phi = pi/4;% Phi estimated for the initial guess     
        Faut = Faut .* exp(-1i*k0*(x_mesh*sin(theta)*cos(phi) + y_mesh*sin(theta)*sin(phi)));
        
        % Evolutionary method. "Global search"       
    elseif Method == 4
       [~, ~, f_mesh]= getFieldLayer(X, Y, Z, Ex, layers(1)); 
       bmi = im2bw(abs(f_mesh)/max(max(abs(f_mesh))),1/sqrt(2));
       x_center = mean(x_mesh(bmi==1)); % Orientation center
       y_center = mean(y_mesh(bmi==1)); % Orientation center
       [TH,PHI,~] = cart2sph(x_center,y_center,layers(1)) ;
       %TH = pi/3; PHI = pi/6;
       k0 = 2*pi/lambda;
       lambda = 0.005;
       I = exp(-1i*k0*(x_mesh*sin(TH)*cos(PHI) + y_mesh*sin(TH)*sin(PHI)));
       Faut = I*1 ;
       
       n = 120*pi;
       k = k0;
       [theta,~,r] = cart2sph(x_mesh,y_mesh,layers(1)) ;
       Er = n*cos(theta)*(1+1./(1j*k*r)).*exp(-1j*k*r)./(2*pi*r.^2) ;
       Ephi = 0;
       Eth = 1j*k*n*sin(theta)*(1+1./(1j*k*r) -1./(k^2*r.^2)).*exp(-1j*k*r)./(4*pi*r) ;
       [Ex_temp,~,~] = sph2cart(Eth, Ephi, Er);
       
       Faut = Ex_temp;
    else
        Faut =  x_mesh *0;    
    end
    
    switch nargin
        case 7
            figure('Name','Initial Guess')
            %figure('Name','Initial Guess','units','normalized','outerposition',[0 0 1 1])
            %Plot initial guess
            subplot(1,2,1); surf(x_mesh,y_mesh, (abs(Faut))); title('Magnitude AUT');   colorbar ; shading interp; view(0,90);
            xlabel('x [mm]'); ylabel('y [mm]');
            subplot(1,2,2); surf(x_mesh,y_mesh, angle(Faut));  title('\angle AUT');   colorbar ; shading interp; view(0,90);
            xlabel('x [mm]'); ylabel('y [mm]');
            if Method == 1
                subplot(1,2,1);  title('Magnitude: AUT - Basic approach');
                subplot(1,2,2);  title('\angle: AUT - Basic approach');
            elseif Method == 2
                subplot(1,2,1);  title(['Magnitude AUT - Nearest Layer z=' num2str(layer_aut) 'mm']);
                subplot(1,2,2);  title(['\angle AUT - Nearest Layer z=' num2str(layer_aut) 'mm']);
            end
    end


end