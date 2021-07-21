%Propagation at 10GHz -- Horn Antenna
clear all
% Simulated Data
Edata = load('Data/EField_z1mm_xy_0.5mm.txt'); % Data of layers from 0mm to 158mm in z

clc;
X=Edata(:,1);
Y=Edata(:,2);
Z=Edata(:,3); %Z = round(Z - 0.8865);

Ex = Edata(:,4)+ 1i*Edata(:,5); % Electrical Field: X component
Ey = Edata(:,6)+ 1i*Edata(:,7); % Electrical Field: Y component
Ez = Edata(:,8)+ 1i*Edata(:,9); % Electrical Field: Z component

%%

f = 60e9; % Frequency of the signal
lambda = 299792458/f; % Wavelength of the signal
layers = [0, 50]; % Distance from the aperture in mm of the layers 1,2
zref = 41.2;
[x_mesh, y_mesh, f_mesh]= getFieldLayer(X, Y, Z, Ex, layers, zref);

Field = f_mesh(:,:,1);
Desired_Field = f_mesh(:,:,2);

dx = 0.0005; % planar x scan step size [m]
dy = 0.0005; % planar y scan step size [m]

E_x = calculatePropagationMatrix(Field, layers, lambda, dx, dy);

%[E_x, x, y] = BackProjection_PlanarNearField_v2(Field,z_0,f,delta_x,delta_y);


      figure; 
      subplot(1,3,1)
        surf(x_mesh, y_mesh, abs(Field)); shading interp; colorbar; view(90,90)
        title(['Input - Max value = ' num2str(max(max(abs(Field))))])
      subplot(1,3,2)
        surf(x_mesh, y_mesh, abs(E_x)); shading interp; colorbar; view(90,90)
        title(['Output - Est. Mag' num2str(max(max(abs(E_x))))])
      subplot(1,3,3)
        surf(x_mesh, y_mesh, abs(Desired_Field)); shading interp; colorbar; view(90,90)
        title(['Desired. Mag. - Max value = ' num2str(max(max(abs(Desired_Field))))])
        
      figure; 
      subplot(1,3,1)
        surf(x_mesh, y_mesh, abs(Field)); shading interp; colorbar; view(180,0)
        title(['Input - Max value = ' num2str(max(max(abs(Field))))])
      subplot(1,3,2)
        surf(x_mesh, y_mesh, abs(E_x)); shading interp; colorbar; view(180,0)
        title(['Output - Est. Mag' num2str(max(max(abs(E_x))))])
      subplot(1,3,3)
        surf(x_mesh, y_mesh, abs(Desired_Field)); shading interp; colorbar; view(180,0)
        title(['Desired. Mag. - Max value = ' num2str(max(max(abs(Desired_Field))))]) 
          
    figure; 
      subplot(1,3,1)
        surf(x_mesh, y_mesh, angle(Field)); shading interp; colorbar; view(90,90)
        title(['Input '])
      subplot(1,3,2)
        surf(x_mesh, y_mesh, angle(E_x)); shading interp; colorbar; view(90,90)
        title(['Output '])
      subplot(1,3,3)
        surf(x_mesh, y_mesh, angle(Desired_Field)); shading interp; colorbar; view(90,90)
        title(['Desired Phase'])
        
%Measured Data