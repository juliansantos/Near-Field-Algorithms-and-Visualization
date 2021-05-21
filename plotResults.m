function plotResults(f_mesh, x_mesh, y_mesh, temp, error, V, i)

    if  V ~=0 
        figure('Name','Simulated and estimated values using PTP algorithm','units','normalized','outerposition',[0 0 1 1])
        plotIterationsIFT(i*1.0, x_mesh, y_mesh, f_mesh, temp)
    end 
    
    figure; plot(error(1,:)); title('Error Layer1 per iterations ' )
    figure; plot(error(2,:)); title('Error Layer2 per iterations ' )

    % Difference in phase layer 1 and 2: Simulated data
        delta_ang_sim_l12 = wrapToPi(angle(f_mesh(:,:,1))-angle(f_mesh(:,:,2))- angle(f_mesh(51,53,1))+angle(f_mesh(52,53,2))); 
    % Difference in magnitude layer 1 and 2: Simulated data    
        delta_mag_sim_l12 = abs(f_mesh(:,:,1)-f_mesh(:,:,2));
    % Difference in phase layer 1 and 2: Estimated data
        delta_ang_est_l12 = wrapToPi(angle(temp(:,:,1))-angle(temp(:,:,2))-angle(temp(52,53,1))+angle(temp(52,53,2)));
    % Difference in magnitud layer 1 and 2: Estimated data
        delta_mag_est_l12 = abs(temp(:,:,1)-temp(:,:,2));
    % 
    
    
    figure('Name','Difference in Magnitude and Phase between layers for Sim. and Est. Values ','units','normalized','outerposition',[0 0 1 1])
    subplot(2,3,1); cla; surf(x_mesh,y_mesh, delta_mag_sim_l12);  title(['Simulation: | \Delta(L1, L2) |  iter:' num2str(i)]);   colorbar ; shading interp;  view(0,90);
    subplot(2,3,4); cla; surf(x_mesh,y_mesh, delta_ang_sim_l12);  title(['Simulation: \Delta (\angle L1, \angle L2)  iter:' num2str(i)]);   colorbar ; shading interp;  view(0,90);
    subplot(2,3,2); cla; surf(x_mesh,y_mesh, delta_mag_est_l12);  title(['Estimation: | \Delta(L1, L2) |  iter: ' num2str(i)]);   colorbar ; shading interp;  view(0,90);
    subplot(2,3,5); cla; surf(x_mesh,y_mesh, delta_ang_est_l12);  title(['Estimation: \Delta(\angle L1, \angle L2)  iter: ' num2str(i)]);  colorbar ; shading interp;  view(0,90);   
    subplot(2,3,3); cla; surf(x_mesh,y_mesh, delta_mag_sim_l12-delta_mag_est_l12);  title(['Sim. | \Delta(L1, L2) | - Est. | \Delta(L1, L2) |' num2str(i)]);   colorbar ; shading interp;  view(0,90);
    subplot(2,3,6); cla; surf(x_mesh,y_mesh, wrapToPi(delta_ang_sim_l12-delta_ang_est_l12));  title(['Sim. \Delta (\angle L1, \angle L2) - Est. \Delta (\angle L1, \angle L2) iter:' num2str(i)]);  colorbar ; shading interp;  view(0,90);  
end 