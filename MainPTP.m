% Author: Julián Santos
% Description: Software written and tested using MatLab2020b, 
%              This code has the following functionalities for a Horn
%              Antenna
%               1. Visualization of the Electric Field, Magnetic Field, Far
%               Field, Power flow. 

%% Preallocation of Variables and Loading Simulated Data
clear all; clc; close all force; 
 
% Loading data from CST 
Edata = load('Data/antBX1mm.txt'); 
%Edata = load('Data/3layerSIM_1mm.txt'); 
    % Elecrical Field Columns:  x [mm], y [mm], z [mm],
    %                           ExRe [V/m], ExIm [V/m], 
    %                           EyRe [V/m], EyIm [V/m], 
    %                           EzRe [V/m], EzIm [V/m].
        
%Hdata = load('Data/HField_z1mm_xy_1mm.txt'); 
    % Magnetic Field Columns:   x[mm], y[mm], z[mm],
    %                           HxRe [A/m], HxIm [A/m], 
    %                           HyRe [A/m], HyIm [A/m], 
    %                           HzRe [A/m], HzIm [A/m].
    
%Pdata= load('Data/PowerFlow_z1mm_xy_1mm.txt'); 
    % Power Flow Columns: x [mm], y [mm], z [mm], 
    %                     Px [V.A/m^2], 
    %                     Py [V.A/m^2], 
    %                     Pz [V.A/m^2].

%FFdata = load('Data/FarField60GHz_Horn.txt'); 
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



% Extraction of parameters
%Simulations
clc;
X=Edata(:,1);
Y=Edata(:,2);
Z=Edata(:,3); %Z = round(Z - 0.8865);
%v = (X == 99 |X == 102 |X == -102 |X == -99);
%X(v)= [];
%Y(v)= [];
%Z(v)=[];


Ex = Edata(:,4)+ 1i*Edata(:,5); % Electrical Field: X component
Ey = Edata(:,6)+ 1i*Edata(:,7); % Electrical Field: Y component
Ez = Edata(:,8)+ 1i*Edata(:,9); % Electrical Field: Z component
%Ex(v)= [];
%Ey(v)= [];
%Ez(v)=[];

%Hx = Hdata(:,4)+ 1i*Hdata(:,5); % Magnetic Field: X component
%Hy = Hdata(:,6)+ 1i*Hdata(:,7); % Magnetic Field: Y component
%Hz = Hdata(:,8)+ 1i*Hdata(:,9); % Magnetic Field: Z component
%Px = Pdata(:,4); % Power: X component
%Py = Pdata(:,5); % Power: X component
%Pz = Pdata(:,6); % Power: X component

NumberEPlanes = length(unique(Edata(:,3))); % Number of E planes = 94
    % The distance between planes is equal to 1mm 
    % The first plane is at z=-0.9 (plane xy)
    % The last plane is at z=92 (plane xy)
NumberHPlanes = length(unique(Edata(:,1))); % Number of H planes = 53
    % The distance between planes is equal to 1mm 
    % The first plane is at x=-26 (plane yz)
    % The last plane is at x=26 (plane yz)

Sim = 1; % Varible to indicate simulated data. 
size_layers = [200, 200; 200, 200];

%% Loading Experimental data and parameter extraction
Edata = load('Data/Data_3layers.txt');
X = Edata(:,1);
Y = Edata(:,2);
zref = 167.5; 
Z = Edata(:,3) - zref; %
Exm = Edata(:,7);
Exm = db2mag(Exm); 
Eym = Exm*0; Ezm = Eym; 
Sim = 0; % Variable to indicate measured data
% Input parameters visualization

    layers = [3,5]; % Distance in milimeters in z direction from the aperture at   
                    % which the layers and fields want to be visualized or obtained.  
    size_layers = [50, 50; 50, 50];  % Dimension in milimeters of the layers 
                                    %to be ploted at the 3D figure
                                    
%clear Edata Hdata Pdata SCdat a; clc;
 
%% Visualization Far Field
    % TO-DO Implement nargin et nargout if neccesary
           % Options: simple, rectangular, slice
           % inside each functions show the functionalities i.e examples of
           % inputs           
           % Explain the argumets use 
    plotFarField(FFdata(:,3),FFdata(:,1),FFdata(:,2));
    
%% Visualization of E field.
    zref = 167.5;
    layers=[20,60];
    %Et = sqrt(abs(Ex).^2+abs(Ey).^2+abs(Ez).^2);
    plotElectricField(X, Y, Z, Ex, Ey, Ez, layers, size_layers);   
    
%% Visualization of H field.  
    layers=[4,27];
    plotMagneticField(X, Y, Z, Hx, Hy, Hz, layers, size_layers); 
    
%% Visualization of Power flow            
    plotPowerFlow(X, Y, Z, Px, Py, Pz, layers, size_layers); 
    
    % Example of power calculated directly from E and H:
        %[Pcx,Pcy,Pcz]=calculatePower(Ex, Ey, Ez, Hx, Hy, Hz);
        %plotPowerFlow(X, Y, Z, Pcx, Pcy, Pcz, layers, size_layers); 

%% Visualization of Impedance
    layers = [10,51];
    impedanceBehaviour(X,Y,Z,Ex, Ey, Ez, Hx, Hy, Hz,layers)   
    
%% PTP algorithm 
% 3 20 50 mm. 
    clc; close all force; 
    f = 10e9; % Frequency of the signal
    lambda = 3e8/f; % Wavelength of the signal
    layers = [20, 40]; % Distance from the aperture in mm of the layers 1,2
                       % for the PTP algorithm 
    cycles = 2000; % Number of iterations for the IFT. 
    zref = 167.5;
                       
    %Please select the type of visualization
        V = 2;
        % Visualization options (Variable V)
            % 0: Show animation convergence phase and magnitude, error best
            % case scenario. If the input data is simulated it will show
            % the best approximation in terms of the error in phase and
            % magnitude, but if it is measured data it will show it just
            % for the first one. Also it shows the delta difference between
            % layers. 
            
            % 1: Show animation just for the convergence of the phase, 
            % also it shows the best case scenarios and deltas estimated
            % and simulated/measured. 
            
            % 2: Show final results without animation: Use this option if
            % you want to work with a remarkable amount of iterations in
            % order to obtain the results faster. 
        
    %Please select the type of initial guess
        I = 0;
        % Initial guess option
            % 0: Not initial guess at Maut -> the algorithm start
            % working with the iterations between layers
            
            % 1: Basic Approach
            
            % 2: Nearest layer above measured/ simulated 
            % to the antenna aperture to propagate it to the first layer
            
            % 3: Hybrid method
            
            % 4: BMI presented by Li Xiang in
            % his master thesis (See article li_xiang_IFT.pdf)
             

    PTP(X, Y, Z, Ex, layers, V, I, cycles, lambda, zref);

    
%% Compare data
clc; close all force;
% Measurement data
Edata = load('Data/Data_3layers.txt');
X = Edata(:,1);
Y = Edata(:,2);
zref = 119; 
Z = Edata(:,3) - zref + 41.2; %
Exm = Edata(:,7);

% Simulated data
Edata = load('Data/3layerSIM_1mm.txt');
X1 = Edata(:,1);
Y1 = Edata(:,2);
zref = 119; 
Z1 = Edata(:,3) - zref + 41.2; %
Exs = Edata(:,4)+ 1i*Edata(:,5);

layers = [3,5];
zref  = 41.2;
[x_mesh, y_mesh, f_mesh1]= getFieldLayer(X, Y, Z, Exm, layers, zref);

[~, ~, f_mesh2]= getFieldLayer(X, Y, Z, abs(Exs), layers, zref);

f_mesh1 = 10.^((f_mesh1)/20);
deltal1 = f_mesh2(:,:,1) - f_mesh1(:,:,1);
deltal2 = f_mesh2(:,:,2) - f_mesh1(:,:,2);
% Compare magnitudes in linear (RAW)
    figure('Name','Linear: Magnitude in layers from Sim. and Mea. Values ','units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,1); cla; surf(x_mesh,y_mesh, f_mesh2(:,:,1));  title('Simulation: | L1 |  ');   colorbar ; shading interp;  view(0,0);
    subplot(2,2,3); cla; surf(x_mesh,y_mesh, f_mesh2(:,:,2));  title('Simulation: | L2 |  ');   colorbar ; shading interp;  view(0,0);
    subplot(2,2,2); cla; surf(x_mesh,y_mesh, f_mesh1(:,:,1));  title('Measurement: | L1 |  ');   colorbar ; shading interp;  view(0,0);
    subplot(2,2,4); cla; surf(x_mesh,y_mesh, f_mesh1(:,:,2));  title('Measurement: | L2 |  ');  colorbar ; shading interp;  view(0,0);  
% Compare magnitudes in dB (RAW)
    figure('Name','dB: Magnitude in layers from Sim. and Mea. Values ','units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,1); cla; surf(x_mesh,y_mesh, mag2db(f_mesh2(:,:,1)));  title('Simulation: | L1 |  ');   colorbar ; shading interp;  view(0,0);
    subplot(2,2,3); cla; surf(x_mesh,y_mesh, mag2db(f_mesh2(:,:,2)));  title('Simulation: | L2 |  ');   colorbar ; shading interp;  view(0,0);
    subplot(2,2,2); cla; surf(x_mesh,y_mesh, mag2db(f_mesh1(:,:,1)));  title('Measurement: | L1 |  ');   colorbar ; shading interp;  view(0,0);
    subplot(2,2,4); cla; surf(x_mesh,y_mesh, mag2db(f_mesh1(:,:,2)));  title('Measurement: | L2 |  ');  colorbar ; shading interp;  view(0,0);  


f_mesh1(:,:,1) = f_mesh1(:,:,1)./max(max(f_mesh1(:,:,1))); 
f_mesh1(:,:,2) = f_mesh1(:,:,2)./max(max(f_mesh1(:,:,2))); 
f_mesh2(:,:,1) = f_mesh2(:,:,1)./max(max(f_mesh2(:,:,1))); 
f_mesh2(:,:,2) = f_mesh2(:,:,2)./max(max(f_mesh2(:,:,2))); 
deltal1 = abs(f_mesh2(:,:,1) - f_mesh1(:,:,1)) ;
deltal2 = abs(f_mesh2(:,:,2) - f_mesh1(:,:,2));

% Compare magnitudes in linear (Normalized)
    figure('Name','Nor: Magnitude in layers from Sim. and Mea. Values ','units','normalized','outerposition',[0 0 1 1])
    subplot(2,3,1); cla; surf(x_mesh,y_mesh, f_mesh2(:,:,1));  title('Simulation: | L1 |  ');   colorbar ; shading interp;  view(0,0);
    subplot(2,3,4); cla; surf(x_mesh,y_mesh, f_mesh2(:,:,2));  title('Simulation: | L2 |  ');   colorbar ; shading interp;  view(0,0);
    subplot(2,3,2); cla; surf(x_mesh,y_mesh, f_mesh1(:,:,1));  title('Measurement: | L1 |  ');   colorbar ; shading interp;  view(0,0);
    subplot(2,3,5); cla; surf(x_mesh,y_mesh, f_mesh1(:,:,2));  title('Measurement: | L2 |  ');  colorbar ; shading interp;  view(0,0);
    subplot(2,3,3); cla; surf(x_mesh,y_mesh, deltal1);  title('Delta | L1 | Sim - Mes  ');   colorbar ; shading interp;  view(0,0);
    subplot(2,3,6); cla; surf(x_mesh,y_mesh, deltal2);  title('Delta | L2 |  Sim - Mes');  colorbar ; shading interp;  view(0,0);  

    
%% Results in phase
    I = 0; V=2; cycles = 20; lambda = 5e-3;
    [x_mesh, y_mesh, f_mesh, temp]=PTP(X, Y, Z, Exs, layers, V, I, cycles, lambda);
    
    [~, ~, f_mesh1, temp2]=PTP(X, Y, Z, db2mag(Exm), layers, V, I, cycles, lambda);
    
    cx = floor(size(x_mesh,1)/2); cy = cx;
    anglel1s = wrapToPi(angle(temp(:,:,2))-angle(temp(cx,cx+1,2)));
    anglel1m = wrapToPi(angle(temp2(:,:,2))-angle(temp2(cx,cx+1,2))); 
    
    subplot(1,3,1); cla; surf(x_mesh,y_mesh, anglel1s);     colorbar ; shading interp;  view(0,90);
    subplot(1,3,2); cla; surf(x_mesh,y_mesh, anglel1m);     colorbar ; shading interp;  view(0,90);
    subplot(1,3,3); cla; surf(x_mesh,y_mesh, wrapToPi(anglel1s - anglel1m));   colorbar ; shading interp;  view(0,90);
    
%% Near to Far Field Transform 
%     clc; close all;
%     %First convert to spheric coordinates to polar coordinates
%     
%     ro = sqrt(abs(X).^2+abs(Y).^2+abs(Z).^2) ;
%     theta = (rad2deg(atan(Y./X)));
%     phi = (rad2deg(atan(sqrt(abs(X).^2+abs(Y).^2)./Z)));
%     
%     E_r = sin(theta).*cos(phi).*Ex + sin(theta).*sin(phi).*Ey +cos(theta).*Ez; 
%     E_theta = cos(theta).*cos(phi).*Ex +cos(theta).*sin(phi).*Ex -sin(theta).*Ez;
%     E_phi = -sin(phi).*Ex + cos(phi).*Ey;
% 
% 
%     Mag = mag2db(sqrt(abs(E_theta).^2 + abs(E_phi).^2));
%     theta(isnan(theta))=0;
%     Input= [(theta),(phi)] ;
%     [iwant, ia, ic] = unique(Input,'rows') ;
%     plotFarField(Mag(ia),iwant(:,1),iwant(:,2));
%     j= 0;
%     %for i=1:length(Mag)
%     %    K(j) = [Mag(i), theta(i), phi(i)];
%     %    j = j+1;
% 
%     %end
%     
    


