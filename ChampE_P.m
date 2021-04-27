% Author: Julián Santos
% Description: Software written and tested using MatLab2020b, 
%              This code has the following functionalities
%               1. Visualization of the Electric Field, Magnetic Field, Far
%               Field, Power flow. 

%% Preallocation Variables and Loading Experimental and Simulated Data
clear all; clc; close all force; 
 
% Loading data from CST 
Edata = load('Data/EField60GHz_Horn.txt'); 
    % Elecrical Field Columns:  x [mm], y [mm], z [mm],
    %                           ExRe [V/m], ExIm [V/m], 
    %                           EyRe [V/m], EyIm [V/m], 
    %                           EzRe [V/m], EzIm [V/m].
        
Hdata = load('Data/HField60GHz_Horn.txt'); 
    % Magnetic Field Columns:   x[mm], y[mm], z[mm],
    %                           HxRe [A/m], HxIm [A/m], 
    %                           HyRe [A/m], HyIm [A/m], 
    %                           HzRe [A/m], HzIm [A/m].
    
Pdata= load('Data/PowerFlow60GHz_Horn.txt'); 
    % Power Flow Columns: x [mm], y [mm], z [mm], 
    %                     Px [V.A/m^2], 
    %                     Py [V.A/m^2], 
    %                     Pz [V.A/m^2].

FFdata = load('Data/FarField60GHz_Horn.txt'); 
    % Far Field Columns: Theta [deg.], Phi [deg.],
    %                    Abs(Dir.) [dBi], Abs(Theta) [dBi], 
    %                    Phase(Theta) [deg.], Abs(Phi)[ dBi], 
    %                    Phase(Phi) [deg.], Ax.Ratio [dB].    
    
SCdata = load('Data/SurfaceCurrent60GHz_Horn.txt'); 
    % Surface Current Columns: x [mm], y [mm], z [mm], 
    %                          KxRe [A/m], KxIm [A/m], 
    %                          KyRe [A/m], KyIm [A/m], 
    %                          KzRe [A/m], KzIm [A/m], 
    %                          Area [mm^2].

%Loading experimental data


%% Extraction of parameters
%Simulations
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

 
%% Visualization Far Field
    % TO-DO Implement nargin et nargout if neccesary
           % Options: simple, rectangular, slice
           % inside each functions show the functionalities i.e examples of
           % inputs           
           % Explain the argumets use
plotFarField(FFdata(:,3),FFdata(:,1),FFdata(:,2));

%% Visualization of E field.
    %TO-DO
    % 3D use nargin to have the option of showing layers
           % Options: simple, rectangular, slice
           % inside each functions show the functionalities i.e examples of
           % inputs
    % Use nargin and nargout at plotting layers at measurement at structure and main plot 
    % input parameter layers
    % obtain simulation data from exact locations
plotElectricField(X, Y, Z, Ex, Ey,Ez);    
%% Visualization of H field.
    % 3D
    % Once developed in E copie et coller ici pour H 
    quiver3(X, Y, Z, abs(Hx), abs(Hy), abs(Hz),'LineWidth',4,'color','r');
    title('Magnetic Field - Magnitud - 60GHz Horn Antenna')
    axis equal
    hold on; plotstructure3D()

%% Plot Patch 2D
% This will be used later for showing the planes of the electrical field 
% as well af the process of estimation of the phases of the system. 
x = [-1 1 1 -1];
y = [-1 -1 1 1];
patch(x,y,'blue')
x1 = [-0.5 0.5 0.5 -0.5];
y1 = [-0.5 -0.5 0.5 0.5];
patch(x1,y1,'black')


%% Show a coordinate 3D system. 
quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[1;0;0],[0;1;0],[0;0;1]); 
% axis('off') % could be of help for representation.


%% parcour Robot 
    %create 3D model 
    %create the 3 different options PNFR CNFR SNFR
    %set referencial points and size of the planes distance betweens plaens
    % [z1, z2, ..., zn] 
    % distance between points = []
    % output: nice 3D with points and robot
    % and file
    


