% here enter the reduced matrix from zero towards chess


% Please see the articles [1]-[14] to see a mathematical demonstration fo the 
% following code. 

function prop_matrix = calculatePropagationMatrix(x_mesh, y_mesh, field, layers, lambda)


% here is a delta lambda.  === Here is calculed the distance between layers
    lambda_mm = lambda*1e3;
    %k = x_mesh.^2+y_mesh.^2;
    %delta_z = sqrt(k+layers(2)^2) - sqrt(k+layers(1)^2);
    %prop_matrix = exp(-1i*2*pi*delta_z/lambda_mm);

   % field =Ex;
% FFT of the 'measurements, 
    fft_field = fft2(field);
% Here it is the integral of the field in two dimensions. 
    ky=2*pi/lambda_mm; kz=ky;

    prop_matrix = ifft2(fft_field*exp(-1i*kz*(layers(2)-layers(1))));
    %prop_matrix = trapz(trapz(value));
% As an example enter the electrical field in x at a specific layer and see
% the behaviour. 
    
end 