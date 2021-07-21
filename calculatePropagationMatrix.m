% here enter the reduced matrix from zero towards chess


% Please see the articles [1]-[14] to see a mathematical demonstration fo the 
% following code. 

function [Propagated_field] = calculatePropagationMatrix(field, layers, lambda0,dx,dy)
    
    k0 = 2*pi/lambda0; 
    [M,N] = size(field);
    MI = 10*M;
    NI = 10*N;
    m = (-MI/2):1:(MI/2-1);
    n = (-NI/2):1:(NI/2-1);

    kx = 2*pi*m/(MI*dx);   % Define kx axis for plane wave spectrum
    ky = 2*pi*n/(NI*dy);   % Define ky axis for plane wave spectrum
    [ky_grid, kx_grid] = meshgrid(ky,kx);
    kz_grid = sqrt(k0^2 - kx_grid.^2 - ky_grid.^2); % Define kz axis for plane wave spectrum
    z_0= (layers(2)-layers(1))*1e-3;
    method = 1;
    
    
    %Method 1 ifft - fft:
        if method == 1
            fx = ifftshift(ifft2(field,MI,NI));
            fx_z0 = fx.*exp(-1i*kz_grid*z_0).*(imag(kz_grid)==0);
            E_x_Zero_Padded = fft2(ifftshift(fx_z0));
            Propagated_field = E_x_Zero_Padded(1:M,1:N);
        elseif method == 2
            
        end
end 