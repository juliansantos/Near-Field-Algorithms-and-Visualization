% here enter the reduced matrix from zero towards chess


% Please see the articles [1]-[14] to see a mathematical demonstration fo the 
% following code. 

function [ifft_fact, factor] = calculatePropagationMatrix(x_mesh, y_mesh, field, layers, lambda)
    
    xref = 0; yref = 0; zref = 15.2 ;
    ko = 2 * pi/lambda;
    %kx = ko *x_mesh*1e-3;
    %ky = ko *y_mesh*1e-3;

    [azimuth1,elevation1,~] = cart2sph(x_mesh,y_mesh,layers(1)*ones(size(x_mesh))+zref);
    kx = ko *sin(elevation1).*cos(azimuth1);
    ky = ko *sin(elevation1).*sin(azimuth1);
    clear azimuth1, elevation1;   
    d = (layers(2)-layers(1))*1e-3;
    factor = exp(-1i.*sqrt(ko.^2 -kx.^2-ky.^2)*d);
   
    method = 1;
    
    %Method 1 ifft - fft:
        if method == 1
            fft_field = ((fft2(field)));
            ifft_fact = (ifft2((fft_field.*factor)));
            
            %figure;
            %surf(abs(fft_field));
            %figure;
    %Method 2 idft - dft: 
        elseif method == 2
            pws = zeros(size(kx));
            % Here implement the DFT of the function 
            for i=1:size(kx,2)
                for j=1:size(ky,1)
                    pws(i,j) = sum(field(i,j)* exp(1i *(kx(i,j)*x_mesh*1e-3 + ky(i,j)*y_mesh*1e-3)),'all');
                end
            end
            % Here implement the IDFT of the function 
            interfact = pws.* factor; 
            ifft_fact = zeros(size(interfact));
            for i=1:size(x_mesh,2)
                for j=1:size(y_mesh,1)
                    ifft_fact(i,j) = sum(interfact(i,j)* exp(-1i *(kx*x_mesh(i,j)*1e-3 + ky*y_mesh(i,j)*1e-3)),'all');
                end
            end
        end
end 