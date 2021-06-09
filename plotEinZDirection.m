% Code to show evolution of Electrical field 

%%
clc; clear all; 
%Edata = load('Data/EField_z1mm_xy_1mm.txt'); 
Edata = load('Data/field2a.txt'); 
%Edata = load('Data/Efield.txt'); 
    % Elecrical Field Columns:  x [mm], y [mm], z [mm],
    %                           ExRe [V/m], ExIm [V/m], 
    %                           EyRe [V/m], EyIm [V/m], 
    %                           EzRe [V/m], EzIm [V/m].

X=Edata(:,1);
Y=Edata(:,2);
Z=Edata(:,3);
Ex = Edata(:,4)+ 1i*Edata(:,5); % Electrical Field: X component
Ey = Edata(:,6)+ 1i*Edata(:,7); % Electrical Field: Y component
Ez = Edata(:,8)+ 1i*Edata(:,9); % Electrical Field: Z component    

Field_3D = reshape(Ey, [length(unique(X)), length(unique(Y)), length(unique(Z))]);
    
%%
close all; clc; 
figure('Name','Efield one point','units','normalized','outerposition',[0 0 1 1]);
x = floor(length(unique(X))/2);  y =x;
%x = 20; y =20;
line_vector = squeeze(Field_3D(x,y,:));
ref_aut = 88;
z_coor = unique(Z)-ref_aut;
subplot(1,2,1);
plot(z_coor,abs(line_vector),'LineWidth',3); title('Magnitude');
xlabel('z[mm]'); grid on;
subplot(1,2,2); 
plot(z_coor,angle(line_vector),'LineWidth',2); title('Phase - Radians');
xlabel('z[mm]');



