%Just propagation at 10GHz -- Horn Antenna

%load data10GHZ % 51 layers: 0...151
f = 10e9; % Frequency of the signal
lambda = 299792458/f; % Wavelength of the signal
layers = [5, 20]; % Distance from the aperture in mm of the layers 1,2
zref = 167.5;
[x_mesh, y_mesh, f_mesh]= getFieldLayer(X, Y, Z, Ex, layers, zref);


Field = f_mesh(:,:,1);
Desired_Field = f_mesh(:,:,2);

dx = 0.001; % planar x scan step size [m]
dy = 0.001; % planar y scan step size [m]

E_x = calculatePropagationMatrix(field_layer1, layers, lambda, dx, dy);

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