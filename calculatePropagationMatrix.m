% here enter the reduced matrix from zero towards chess


% Please see the articles [1]-[14] to see a mathematical demonstration fo the 
% following code. 

function [prop_matrix, factor] = calculatePropagationMatrix(x_mesh, y_mesh, field, layers, lambda)

    %k = x_mesh.^2+y_mesh.^2;
    %delta_z = sqrt(k+layers(2)^2) - sqrt(k+layers(1)^2);
    %prop_matrix = exp(-1i*2*pi*delta_z/lambda_mm);

    % FFT of the 'measurements/ simulations' 
        fft_field = fft2(field);
    % Here it is the integral of the field in two dimensions. 
    

        ko = 2 * pi/lambda;
        %kx = ko *x_mesh*1e-3;
        %ky = ko *y_mesh*1e-3;
        
        [azimuth1,elevation1,r1] = cart2sph(x_mesh,y_mesh,layers(1)*ones(size(x_mesh)));
        %[azimuth2,elevation2,r2] = cart2sph(x_mesh,y_mesh,layers(2)*ones(size(x_mesh)));
        kx = ko *sin(elevation1).*cos(azimuth1);
        ky = ko *sin(elevation1).*sin(azimuth1);
        %clear r1, azimuth1, elevation1; 
        
        
        d = (layers(2)-layers(1))*1e-3;
        factor = sqrt(ko.^2 -kx.^2-ky.^2)*d;
        prop_matrix = ifft2(fft_field.*exp(-1i.*factor));
    
    
end 