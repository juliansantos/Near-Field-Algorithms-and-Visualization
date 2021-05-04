% Author: Julián Santos
% Description: Software written and tested using MatLab2020b, 
%              This code has the following functionalities for a Horn
%              Antenna
%               1. Visualization of the Electric Field, Magnetic Field, Far
%               Field, Power flow. 

%% Preallocation of Variables and Loading Experimental and Simulated Data
clear all; clc; close all force; 
 
% Loading data from CST 
Edata = load('Data/EField_z1mm_xy_1mm.txt'); 
    % Elecrical Field Columns:  x [mm], y [mm], z [mm],
    %                           ExRe [V/m], ExIm [V/m], 
    %                           EyRe [V/m], EyIm [V/m], 
    %                           EzRe [V/m], EzIm [V/m].
        
Hdata = load('Data/HField_z1mm_xy_1mm.txt'); 
    % Magnetic Field Columns:   x[mm], y[mm], z[mm],
    %                           HxRe [A/m], HxIm [A/m], 
    %                           HyRe [A/m], HyIm [A/m], 
    %                           HzRe [A/m], HzIm [A/m].
    
Pdata= load('Data/PowerFlow_z1mm_xy_1mm.txt'); 
    % Power Flow Columns: x [mm], y [mm], z [mm], 
    %                     Px [V.A/m^2], 
    %                     Py [V.A/m^2], 
    %                     Pz [V.A/m^2].

FFdata = load('Data/FarField60GHz_Horn.txt'); 
    % Far Field Columns: Theta [deg.], Phi [deg.],
    %                    Abs(Dir.) [dBi], Abs(Theta) [dBi], 
    %                    Phase(Theta) [deg.], Abs(Phi)[ dBi], 
    %                    Phase(Phi) [deg.], Ax.Ratio [dB].    
    
%SCdata = load('Data/SurfaceCurrent60GHz_Horn.txt'); 
    % Surface Current Columns: x [mm], y [mm], z [mm], 
    %                          KxRe [A/m], KxIm [A/m], 
    %                          KyRe [A/m], KyIm [A/m], 
    %                          KzRe [A/m], KzIm [A/m], 
    %                          Area [mm^2].

%Loading experimental data


%% Extraction of parameters
%Simulations
clc;
X=Edata(:,1);
Y=Edata(:,2);
Z=Edata(:,3);
Ex = Edata(:,4)+ 1i*Edata(:,5); % Electrical Field: X component
Ey = Edata(:,6)+ 1i*Edata(:,7); % Electrical Field: Y component
Ez = Edata(:,8)+ 1i*Edata(:,9); % Electrical Field: Z component
Hx = Hdata(:,4)+ 1i*Hdata(:,5); % Magnetic Field: X component
Hy = Hdata(:,6)+ 1i*Hdata(:,7); % Magnetic Field: Y component
Hz = Hdata(:,8)+ 1i*Hdata(:,9); % Magnetic Field: Z component
Px = Pdata(:,4); % Power: X component
Py = Pdata(:,5); % Power: X component
Pz = Pdata(:,6); % Power: X component

NumberEPlanes = length(unique(Edata(:,3))); % Number of E planes = 94
    % The distance between planes is equal to 1mm 
    % The first plane is at z=-0.9 (plane xy)
    % The last plane is at z=92 (plane xy)
NumberHPlanes = length(unique(Edata(:,1))); % Number of H planes = 53
    % The distance between planes is equal to 1mm 
    % The first plane is at x=-26 (plane yz)
    % The last plane is at x=26 (plane yz)

% Input parameters visualization

    layers = [-14,30]; % Distance in milimeters in z direction from the aperture at   
                    % which the layers and fields want to be visualized or obtained.  
    size_layers = [50, 50; 50, 50];  % Dimension in milimeters of the layers 
                                    %to be ploted at the 3D figure
                                    
clear Edata Hdata Pdata SCdata; clc;
 
%% Visualization Far Field
    % TO-DO Implement nargin et nargout if neccesary
           % Options: simple, rectangular, slice
           % inside each functions show the functionalities i.e examples of
           % inputs           
           % Explain the argumets use
    plotFarField(FFdata(:,3),FFdata(:,1),FFdata(:,2));
                                    
%% Visualization of E field.         
    layers=[15,38.6];
    plotElectricField(X, Y, Z, Ex, Ey,Ez, layers, size_layers);   
    
%% Visualization of H field.           
    plotMagneticField(X, Y, Z, Hx, Hy, Hz, layers, size_layers); 
    
%% Visualization of Power flow            
    plotPowerFlow(X, Y, Z, Px, Py, Pz, layers, size_layers); 
    
    % Example of power calculated directed from E and H:
        %[Pcx,Pcy,Pcz]=calculatePower(Ex, Ey, Ez, Hx, Hy, Hz);
        %plotPowerFlow(X, Y, Z, Pcx, Pcy, Pcz, layers, size_layers); 
        
%% Visualization of Impedance
    layers = [1,51];
    impedanceBehaviour(X,Y,Z,Ex, Ey, Ez, Hx, Hy, Hz,layers)   
%% PTP algorithm 
    clc
    f = 60e9; % Frequency of the signal
    lambda = 3e8/f; % Wavelength of the signal
    layers = [5, 15]; % Distance in milimeters in z direction 
                     % from the aperture for the PTP algorithm 
    
    %Getting the fields to work with
    [x_mesh, y_mesh, f_mesh]= getFieldLayer(X, Y, Z, Ex, layers);
    %Fields --> M<P   (M=magnitude) (P=Phase)
        %Getting the magnitude of the fields (Input of the algorithm) 
            M=abs(f_mesh);
        %Getting the phases of the fields (for testing the algorithm)
            P=angle(f_mesh);
          
            
     %[x_mesh, y_mesh, g_mesh]= getFieldLayer(X, Y, Z, Hx, layers); 
     %       M1=abs(g_mesh);
     %       P1=angle(g_mesh);
    
    
      %Initial guess enter to function calculate Propagation Matrix and 
      	Maut = initialguess(x_mesh, y_mesh);
      
      cycles = 100;
      error = zeros(1,cycles);
      

      %AUT to Layer 1  
          %Propagating initial guess to first layer
            field_layer1 = calculatePropagationMatrix(x_mesh, y_mesh, Maut, [0 layers(1)], lambda);
          %Replace amplitudes estimated by measurements/simulated
            field_layer1 = M(:,:,1).*exp(1i*angle(field_layer1));
            
     for i =1:cycles           
      %Layer 1 to Layer 2 
          %Propagating layer 1 to layer 2
            field_layer2 = calculatePropagationMatrix(x_mesh, y_mesh, field_layer1, layers, lambda);
          %Replace amplitudes estimated by measurements/simulated
            field_layer2 = M(:,:,2).*exp(1i*angle(field_layer2));
        
      %Layer 2 to Layer 1 
          %Propagating initial guess to first layer
            field_layer1 = calculatePropagationMatrix(x_mesh, y_mesh, field_layer2, flip(layers), lambda);
            
      % check criteria of error field_layer1 vs measurements see behaviour
      % in time -> MSE 
      
            error(i) = sum((abs(field_layer1) - M(:,:,1))./M(:,:,1), 'all');
      
          %Replace amplitudes estimated by measurements/simulated
              
         field_layer1 = M(:,:,1).*exp(1i*angle(field_layer1));   
    
     end 
        close all
        subplot(2,2,1); surf(x_mesh,y_mesh, mag2db(abs(field_layer1))); view(0,90); title('Layer 1 Mag');   colorbar ; shading interp; 
        subplot(2,2,2); surf(x_mesh,y_mesh, angle(field_layer1)); view(0,90); title('Layer 1 Phas');   colorbar ; shading interp; 
        subplot(2,2,3); surf(x_mesh,y_mesh, mag2db(abs(field_layer2))); view(0,90); title('Layer 2 Mag');   colorbar ; shading interp; 
        subplot(2,2,4); surf(x_mesh,y_mesh, angle(field_layer2)); view(0,90); title('Layer 2 Phas');   colorbar ; shading interp; 

        figure; plot(error);
        

        
        %close all
        %subplot(2,2,1); surf(x_mesh,y_mesh, angle((f_mesh(:,:,1)))); view(0,90); title('Mesh 0');   colorbar ; 
        %subplot(2,2,2); surf(x_mesh,y_mesh, angle(f_mesh(:,:,2))); view(0,90); title('Mesh 1');   colorbar ; 
        %subplot(2,2,3); surf(x_mesh,y_mesh, angle(prop_matrix)); view(0,90); title('Estimated mesh 1');    colorbar ;
        
        %See the behaviour of the phase  
        %close all 
        %a = abs(rad2deg(P(:,:,1) - P(:,:,2)))-72;
        %subplot(2,4,1); surf(x_mesh,y_mesh, angle(P(:,:,1)*exp(-1i*2*pi*0))); view(0,90);    colorbar ;
        %subplot(2,4,2); surf(x_mesh,y_mesh, angle(P(:,:,1)*exp(-1i*2*pi*0.1))); view(0,90);    colorbar ;
        %subplot(2,4,3); surf(x_mesh,y_mesh, angle(P(:,:,1)*exp(-1i*2*pi*0.3))); view(0,90);    colorbar ;
        %subplot(2,4,4); surf(x_mesh,y_mesh, angle(P(:,:,1)*exp(-1i*2*pi*2))); view(0,90);    colorbar ;
        %subplot(2,4,5); surf(x_mesh,y_mesh, angle(P(:,:,1)*exp(-1i*2*pi*3))); view(0,90);    colorbar ;
        %subplot(2,4,6); surf(x_mesh,y_mesh, angle(P(:,:,1)*exp(-1i*2*pi*3.1))); view(0,90);    colorbar ;
        %subplot(2,4,7); surf(x_mesh,y_mesh, angle(P(:,:,1)*exp(-1i*2*pi*4))); view(0,90);    colorbar ;
        %subplot(2,4,8); surf(x_mesh,y_mesh, angle(P(:,:,1)*exp(-1i*2*pi*5))); view(0,90);    colorbar ;

    
            
        % (angle(exp([1:10]*1e-3*1i*2*pi/lambda)))    
% Obtain layers
% Clear variables not interested in to save memory



    


