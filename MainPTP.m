% Author: Julián Santos
% Description: Software written and tested using MatLab2020b, 
%              This code has the following functionalities for a Horn
%              Antenna
%               1. Visualization of the Electric Field, Magnetic Field, Far
%               Field, Power flow. 

%% Preallocation Variables and Loading Experimental and Simulated Data
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

clear Edata Hdata Pdata SCdata;
 
%% Visualization Far Field
    % TO-DO Implement nargin et nargout if neccesary
           % Options: simple, rectangular, slice
           % inside each functions show the functionalities i.e examples of
           % inputs           
           % Explain the argumets use
    plotFarField(FFdata(:,3),FFdata(:,1),FFdata(:,2));

%% Visualization of E field.
    clc;
    layers = [15,30]; % Distance in milimeters in z direction from the aperture at   
                    % which thelayers and fields want to be visualized. 
    size_layers = [50, 50; 50, 50];  % Dimension in milimeters of the layers 
                                    %to be ploted at the 3D figure            
    plotElectricField(X, Y, Z, Ex, Ey,Ez, layers, size_layers);    
%% Visualization of H field.
    clc;
    layers = [1,15.8]; % Distance in milimeters in z direction from the aperture at   
                    % which thelayers and fields want to be visualized. 
    size_layers = [50, 50; 50, 50];  % Dimension in milimeters of the layers 
                                    %to be ploted at the 3D figure            
    plotMagneticField(X, Y, Z, Hx, Hy, Hz, layers, size_layers); 
    
%% Visualization of Power flow 
    
    clc;
    layers = [1,15.8]; % Distance in milimeters in z direction from the aperture at   
                    % which thelayers and fields want to be visualized. 
    size_layers = [50, 50; 50, 50];  % Dimension in milimeters of the layers 
                                    %to be ploted at the 3D figure            
    plotPowerFlow(X, Y, Z, Px, Py, Pz, layers, size_layers); 
%% PTP algorithm 

% obtain layers
% clear variables not interested in to save memory



    


