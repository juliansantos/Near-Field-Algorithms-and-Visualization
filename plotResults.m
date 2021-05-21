function plotResults(f_mesh,x_mesh, y_mesh,temp1,temp2,error,V,i)
        close all
        figure('units','normalized','outerposition',[0 0 1 1])
        plotIterationsIFT(i*1.0, x_mesh, y_mesh, f_mesh(:,:,1), f_mesh(:,:,2), temp2, temp1)
        figure; plot(error(1,:)); title('Error Layer1 per iterations ' )
        figure; plot(error(2,:)); title('Error Layer2 per iterations ' )
        
        figure;
        subplot(2,3,1); cla; surf(x_mesh,y_mesh, mag2db(abs(f_mesh(:,:,1)-f_mesh(:,:,2))));  title(['|Layer 1| -- Delta Layer 1 Mag: ' num2str(i)]);   colorbar ; shading interp;  view(0,90);
        subplot(2,3,4); cla; surf(x_mesh,y_mesh, angle(f_mesh(:,:,1))-angle(f_mesh(:,:,2)));  title(['\angle Layer 1 -- Estimated Values iter:' num2str(i)]);   colorbar ; shading interp;  view(0,90);
        subplot(2,3,2); cla; surf(x_mesh,y_mesh, mag2db(abs(temp1-temp2)));  title(['|Layer 2| -- Delta Layer 2 Mag: ' num2str(i)]);   colorbar ; shading interp;  view(0,90);
        subplot(2,3,5); cla; surf(x_mesh,y_mesh, angle(temp1)-angle(temp2));  title(['\angle Layer 2 -- Estimated Values iter:' num2str(i)]);  colorbar ; shading interp;  view(0,90);   
        subplot(2,3,3); cla; surf(x_mesh,y_mesh, mag2db(abs(f_mesh(:,:,1)-f_mesh(:,:,2)))-mag2db(abs(temp1-temp2)));  title(['|Layer 2| -- Delta Layer 2 Mag: ' num2str(i)]);   colorbar ; shading interp;  view(0,90);
        subplot(2,3,6); cla; surf(x_mesh,y_mesh, angle(f_mesh(:,:,1))-angle(f_mesh(:,:,2)) - angle(temp1)+angle(temp2));  title(['\angle Layer 2 -- Estimated Values iter:' num2str(i)]);  colorbar ; shading interp;  view(0,90);  
end 