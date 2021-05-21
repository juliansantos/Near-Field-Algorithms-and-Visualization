function plotIterationsIFT(i, x_mesh, y_mesh, M_layer1, M_layer2, ME_layer2, ME_layer1)
    % i: iteration number i
    % M_layer1 : magnitude fields layer 1 -> simulated/measured
    % M_layer2: magnitude fileds layer 2 -> simulated/measured
    % factor: this is the factor that modifies how change the fields from
    % layer one to layer 2
    % ME_layer1: magnitude of the estimated field at layer one from layer 2
    
    if i==1 || ~isinteger( i ) 
        %figure('units','normalized','outerposition',[0 0 1 1])
        subplot(2,4,1); surf(x_mesh,y_mesh, mag2db(abs(M_layer1))); title('|Layer 1| -- Sim. Values');   colorbar ; shading interp; view(0,90);
        subplot(2,4,5); surf(x_mesh,y_mesh, angle(M_layer1));  title('\angle Layer 1-- Sim. Values');   colorbar ; shading interp; view(0,90);
        subplot(2,4,2); surf(x_mesh,y_mesh, mag2db(abs(M_layer2)));  title('|Layer 2| -- Sim. Values');   colorbar ; shading interp; view(0,90);
        subplot(2,4,6); surf(x_mesh,y_mesh, angle(M_layer2)-angle(M_layer2(52,52)));  title('\angle Layer 2 -- Sim. Values');   colorbar ; shading interp;  view(0,90);
    end
    
    
    subplot(2,4,3); cla; surf(x_mesh,y_mesh, mag2db(abs(ME_layer1)));  title(['|Layer 1| -- Estimated Values iter: ' num2str(i)]);   colorbar ; shading interp;  view(0,90);
    subplot(2,4,7); cla; surf(x_mesh,y_mesh, angle((ME_layer1)));  title(['\angle Layer 1 -- Estimated Values iter:' num2str(i)]);   colorbar ; shading interp;  view(0,90);
    subplot(2,4,4); cla; surf(x_mesh,y_mesh, mag2db(abs(ME_layer2)));  title(['|Layer 2| -- Estimated Values iter: ' num2str(i)]);   colorbar ; shading interp;  view(0,90);
    subplot(2,4,8); cla; surf(x_mesh,y_mesh, angle((ME_layer2)));  title(['\angle Layer 2 -- Estimated Values iter:' num2str(i)]);  colorbar ; shading interp;  view(0,90);   
    pause(1e-3)
    
    

end 