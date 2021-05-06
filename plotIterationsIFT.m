function plotIterationIFT(i, x_mesh, y_mesh, M_layer1, M_layer2, factor, ME_layer1)
    % i: iteration number i
    % M_layer1 : magnitude fields layer 1 -> simulated/measured
    % M_layer2: magnitude fileds layer 2 -> simulated/measured
    % factor: this is the factor that modifies how change the fields from
    % layer one to layer 2
    % ME_layer1: magnitude of the estimated field at layer one from layer 2
    
    subplot(2,4,1); surf(x_mesh,y_mesh, mag2db(abs(M_layer1))); title('|Layer 1| -- Real Values');   colorbar ; shading interp; view(0,90);
    subplot(2,4,5); surf(x_mesh,y_mesh, angle(M_layer1));  title('< Layer 1-- Real Values');   colorbar ; shading interp; view(0,90);
    subplot(2,4,2); surf(x_mesh,y_mesh, mag2db(abs(M_layer2)));  title('|Layer 2| -- Real Values');   colorbar ; shading interp; view(0,90);
    subplot(2,4,6); surf(x_mesh,y_mesh, angle(M_layer2));  title('< Layer 2 -- Real Values');   colorbar ; shading interp;  view(0,90);
    subplot(2,4,3); surf(x_mesh,y_mesh, mag2db(abs(ME_layer1)));  title('|Layer 1| -- Estimated Values');   colorbar ; shading interp;  view(0,90);
    subplot(2,4,7); surf(x_mesh,y_mesh, angle(ME_layer1));  title('< Layer 1 -- Estimated Values');   colorbar ; shading interp;  view(0,90);
    subplot(2,4,4); surf(x_mesh,y_mesh, mag2db(abs(factor)));  title('|factor| -- z2 -> z1');   colorbar ; shading interp;  view(0,90);
    subplot(2,4,8); surf(x_mesh,y_mesh, angle(factor));  title('< factor 1 -- z2 -> z1');   colorbar ; shading interp;  view(0,90);   
    pause(0.01)
    
    

end 