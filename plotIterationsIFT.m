function plotIterationsIFT(i, x_mesh, y_mesh, M_layer, ME_layer)
    % i: iteration number i
    % M_layer1 : magnitude fields layer 1 -> simulated/measured
    % M_layer2: magnitude fileds layer 2 -> simulated/measured
    % factor: this is the factor that modifies how change the fields from
    % layer one to layer 2
    % ME_layer1: magnitude of the estimated field at layer one from layer 2
    % Center coordinate to normalize
    cx = floor(size(x_mesh,1)/2); 
    
    if i==1 || ~isinteger( i ) 
        %figure('units','normalized','outerposition',[0 0 1 1])
        subplot(2,4,1); surf(x_mesh,y_mesh, (abs(M_layer(:,:,1)))); title(['|L1| -- Sim. Value - Max value = ' num2str(max(max(abs(M_layer(:,:,1)))))]);   colorbar ; shading interp; view(0,90);
        subplot(2,4,5); surf(x_mesh,y_mesh, wrapToPi(wrapToPi(angle(M_layer(:,:,1)) - angle(M_layer(cx,cx,1) ))));  title('\angle L1-- Sim. Values');   colorbar ; shading interp; view(0,90);
        subplot(2,4,2); surf(x_mesh,y_mesh, (abs(M_layer(:,:,2))));  title(['|L2| -- Sim. Value - Max value = ' num2str(max(max(abs(M_layer(:,:,2)))))]);   colorbar ; shading interp; view(0,90);
        subplot(2,4,6); surf(x_mesh,y_mesh, wrapToPi(wrapToPi(angle(M_layer(:,:,2)) - angle(M_layer(cx,cx,2)))));  title('\angle L2 -- Sim. Values');   colorbar ; shading interp;  view(0,90);
    end
    
    
    subplot(2,4,3); cla; surf(x_mesh,y_mesh, (abs(ME_layer(:,:,1))));  title(['|L1| -- Est. Value - Max value = ' num2str(max(max(abs(ME_layer(:,:,1))))) ' iter: ' num2str(i)]);   colorbar ; shading interp;  view(0,90);
    subplot(2,4,7); cla; surf(x_mesh,y_mesh, wrapToPi(wrapToPi(angle(ME_layer(:,:,1))-angle(ME_layer(cx,cx,1))) ) );  title(['\angle L1 -- Estimated Values iter:' num2str(i)]);   colorbar ; shading interp;  view(0,90);
    subplot(2,4,4); cla; surf(x_mesh,y_mesh, (abs(ME_layer(:,:,2))));  title(['|L2| -- Est. Value - Max value = ' num2str(max(max(abs(ME_layer(:,:,2)))))  ' iter: ' num2str(i)]);   colorbar ; shading interp;  view(0,90);
    subplot(2,4,8); cla; surf(x_mesh,y_mesh, wrapToPi(wrapToPi(angle((ME_layer(:,:,2))) -angle(ME_layer(cx,cx,2))))  );  title(['\angle L2 -- Estimated Values iter:' num2str(i)]);  colorbar ; shading interp;  view(0,90);   
    pause(1e-3)
    
    

end 