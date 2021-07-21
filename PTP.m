function [x_mesh, y_mesh, f_mesh, temp]=PTP(X, Y, Z, Ex, layers, V, I, cycles, lambda, zref, dx, dy) 

        %Getting the fields to work with (Layer 1 and 2)
    [x_mesh, y_mesh, f_mesh]= getFieldLayer(X, Y, Z, Ex, layers, zref);
    
    %Magnitude and phase of the Fields --> (M=magnitude) (P=Phase)
        %Getting the magnitude of the fields (Input of the algorithm) 
            M=abs(f_mesh);
        %Getting the phases of the fields (for testing the algorithm)
            P=angle(f_mesh);
          
    %Initial guess enter to function calculate Propagation Matrix and 
     [Maut, layer_aut] = initialguess(X, Y, Z, Ex, I, layers,'plot', zref);
            
      %Preallocating and starting variables for iterations 
        
        error = zeros(2,cycles); % Error in Magnitude on each layer per iteration.
        
        errorp = zeros(2,cycles); % Error in Phase in each layer per iteration.
           % This is just used in the case of simulation for the sake of
           % demonstrating the convergence to the true phase 
        
        % best case variables with respect to magnitude error
        best_case = zeros([size(Maut),2]);    
        error_best_case = [inf, inf]; 
        iter_best_case = [1, 1];
        
        % best case variables with respect to phase error "just for the case of simulated data"
            % This is aimed to show the convergence behavior of the 'real'
            % phase with the control parameter of the magnitude of the
            % electric field.           
        best_case_p = zeros([size(Maut),2]);    
        error_best_case_p = [inf, inf]; 
        iter_best_case_p = [1, 1];
        
       %[~,~,Maut] = getFieldLayer(X, Y, Z, Ex, layers(1));  Maut = abs(Maut);
       
      %AUT to Layer 1  
          %Propagating initial guess to first layer
            if I~=0 
                field_layer1 = calculatePropagationMatrix(Maut, [layer_aut layers(1)], lambda,dx,dy); 
            else
                field_layer1 = 0;
            end
          %Replace amplitudes estimated by measurements/simulated
            field_layer1 = M(:,:,1).*exp(1i*angle(field_layer1));            
        
     if V==0, figure('Name','Simulated and estimated values using PTP algorithm (Animation)','units','normalized','outerposition',[0 0 1 1]),end
     if V==1, figure('Name','Estimated phase values using PTP algorithm (Animation)','units','normalized','outerposition',[0 0 1 1]),end
     
     i=1;% variable for the iterations
     flag = 1; % variable for stop criteria
     threshold = 1e-3;
     while  flag~=0           
      %Layer 1 to Layer 2 
          %Propagating layer 1 to layer 2
            field_layer2 = calculatePropagationMatrix(field_layer1, layers, lambda,dx,dy); temp(:,:,2) = field_layer2;
          %Replace amplitudes estimated by measurements/simulated
            field_layer2 = M(:,:,2).*exp(1i*angle(field_layer2));
        
      %Layer 2 to Layer 1 
          %Propagating initial guess to first layer
             field_layer1 = calculatePropagationMatrix(field_layer2, flip(layers), lambda,dx,dy); temp(:,:,1) = field_layer1;
          %Replace amplitudes estimated by measurements/simulated 
             field_layer1 = M(:,:,1).*exp(1i*angle(field_layer1));
           
          %Animation of estimated values using PTP algorithm
             if V==0, plotIterationsIFT(i, x_mesh, y_mesh, f_mesh, temp), end
          %Animation only for phase convergence
             if V==1
                 subplot(1,2,1); cla; surf(x_mesh,y_mesh, wrapToPi(P(:,:,1)-P(52,53,1)));  title(['\angle Layer 1 -- Simulated Values iter:' num2str(i)]);   colorbar ; shading interp;  view(0,90);
                 subplot(1,2,2); cla; surf(x_mesh,y_mesh, wrapToPi(angle(field_layer1)-angle(field_layer1(52,53))));  title(['\angle Layer 1 -- Estimated Values iter:' num2str(i)]); colorbar ; shading interp;  view(0,90), pause(1e-3);
             end
             
      %Error in magnitude0 and phase 
           % error(1,i) = sum((abs(temp(:,:,1)) - M(:,:,1)).^2, 'all')/sum(M(:,:,1).^2,'all');
           % error(2,i) = sum((abs(temp(:,:,2)) - M(:,:,2)).^2, 'all')/sum(M(:,:,2).^2,'all');
           % errorp(1,i) = sum(wrapToPi(angle(temp(:,:,1)) - P(:,:,1)).^2, 'all')/sum(P(:,:,1).^2,'all');
           % errorp(2,i) = sum(wrapToPi(angle(temp(:,:,2)) - P(:,:,2)).^2, 'all')/sum(P(:,:,2).^2,'all');
      
            error(1,i) = sum((abs(temp(:,:,1)) - M(:,:,1)).^2, 'all');
            error(2,i) = sum((abs(temp(:,:,2)) - M(:,:,2)).^2, 'all');
            errorp(1,i) = sum(wrapToPi(angle(temp(:,:,1)) - P(:,:,1)).^2, 'all');
            errorp(2,i) = sum(wrapToPi(angle(temp(:,:,2)) - P(:,:,2)).^2, 'all');
            
      %Best case scenario in magnitude error layer 1       
            if error_best_case(1) > error(1,i) 
                best_case(:,:,1) = temp(:,:,1);
                error_best_case(1) = error(1,i);
                iter_best_case(1) = i;
            end
       %Best case scenario in magnitude error layer 2       
            if error_best_case(2) > error(2,i) 
                best_case(:,:,2) = temp(:,:,2);
                error_best_case(2) = error(2,i);
                iter_best_case(2) = i;
            end      
            
       %Best case scenario in phase error layer 1       
            if error_best_case_p(1) > errorp(1,i) 
                best_case_p(:,:,1) = temp(:,:,1);
                error_best_case_p(1) = errorp(1,i);
                iter_best_case_p(1) = i;
            end
       %Best case scenario in phase error layer 2       
            if error_best_case_p(2) > errorp(2,i) 
                best_case_p(:,:,2) = temp(:,:,2);
                error_best_case_p(2) = errorp(2,i);
                iter_best_case_p(2) = i;
            end    
            
      %Stop criteria based on number of iterations or threshold error in magnitude      
            if error(1,i)<threshold || i>=cycles 
               flag = 0;
            end
            i = i+1;
     end 
        % Plot results
        plotResults(f_mesh,x_mesh, y_mesh, temp, error, errorp, best_case, best_case_p, V, cycles, iter_best_case, iter_best_case_p);
        
end 