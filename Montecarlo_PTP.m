%Description: Montecarlo Simulations -> Propagation Function
    clc;
    lambda = 3e8/10e9;
    layers = [3, 7]; % milimeters. 
    maxval = 50; % Max limit of planes
    
    ko = 2 * pi/lambda;
    
    V = (-maxval:maxval);
    x_mesh = repmat(V,length(V),1);
    y_mesh = x_mesh';
    Nmeasures = length(V)^2; 
    [azimuth1,elevation1,~] = cart2sph(x_mesh,y_mesh,layers(2)*ones(size(x_mesh)));
    [azimuth2,elevation2,~] = cart2sph(x_mesh,y_mesh,layers(1)*ones(size(x_mesh)));
    kx1 = ko *sin(elevation1).*cos(azimuth1);
    ky1 = ko *sin(elevation1).*sin(azimuth1); 
    kx2 = ko *sin(elevation2).*cos(azimuth2);
    ky2 = ko *sin(elevation2).*sin(azimuth2);  
    d = (layers(2)-layers(1))*1e-3;
    factor1 = exp(-1i.*sqrt(ko.^2 -kx1.^2-ky1.^2)*d);
    factor2 = exp(1i.*sqrt(ko.^2 -kx2.^2-ky2.^2)*d);
    %
    maxval_field1 = sqrt(3000);
    maxval_field2 = sqrt(2500);
    
    iter2 = 100;
    iter = 100;
    % Preallocation;
    change1 = zeros(iter,iter2); change2 = change1; 
    error1 = change1; error2 = change1; 
    for j = 1:iter2
        field1 = maxval_field1*rand(size(x_mesh))+ 1i*maxval_field1*rand(size(x_mesh));
        field2 = maxval_field2*rand(size(x_mesh))+ 1i*maxval_field2*rand(size(x_mesh));

        M1 = abs(field1);
        M2 = abs(field2);

        for i=1:iter
            % Layer 1 to layer 2
            fft_field = ((fft2(field1)));
            ifft_fact = (ifft2((fft_field.*factor1)));
            change1(i,j) = sum((M1-abs(ifft_fact))^2, 'all')/ Nmeasures;
            error1(i,j) = sum((M2-abs(ifft_fact))^2, 'all')/ Nmeasures;

            field2 = M2.*exp(1i*angle(ifft_fact)); % Replace Magnitudes

            fft_field = ((fft2(field2)));
            ifft_fact = (ifft2((fft_field.*factor2)));
            change2(i,j) = sum((M2-abs(ifft_fact))^2, 'all')/ Nmeasures;
            error2(i,j) = sum((M1-abs(ifft_fact))^2, 'all')/ Nmeasures;

            field1  = M1.*exp(1i*angle(ifft_fact));
        end
    end

    %% plot results
    figure;
    for j=1:iter2
       plot(100*change1(:,j)/maxval_field1^2) ;
       set(gca,'FontSize',18)
       title('% Av. change Layer1', 'FontSize', 20)
       xlabel('Iterations', 'FontSize', 24)
       ylabel('Change in magnitude % ', 'FontSize', 24)
       hold on;
    end
    
    
    figure;
    for j=1:iter2
       plot(100*change2(:,j)/maxval_field2^2) ;
       set(gca,'FontSize',18)
       title('% Av. change Layer2', 'FontSize', 20)
       xlabel('Iterations', 'FontSize', 24)
       ylabel('Change in magnitude % ', 'FontSize', 24)
       hold on;
    end
    
