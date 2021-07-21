function plotResults(f_mesh, x_mesh, y_mesh, temp, error, errorp, best_case, best_case_p, V, i, iter_best_case, iter_best_case_p)

    if  V ~=0 
        figure('Name','Simulated and estimated values using PTP algorithm (Last Iteration)','units','normalized','outerposition',[0 0 1 1])
        plotIterationsIFT(i*1.0, x_mesh, y_mesh, f_mesh, temp)
    end 
    
    figure; plot(error(1,:)); hold on;
    plot(error(2,:)); title('Error in Magnitud per iterations' )
    %plot(error(1,:)-error(2,:)); 
    xline(iter_best_case(1),'--r');
    xline(iter_best_case(2),'--b');
    legend('Layer1', 'Layer2', 'Minimum error Layer 1', 'Minimum error Layer 2'); grid on;
    
    %figure; plot(errorp(1,:)); hold on;
    %plot(errorp(2,:)); title('Error in Phase per iterations' )
    %plot(errorp(2,:)-error(1,:)); 
    %xline(iter_best_case_p(1),'--r');
    %xline(iter_best_case_p(2),'--b');
    %legend('Layer1', 'Layer2','Difference', 'Minimum error Layer 1', 'Minimum error Layer 2'); grid on;
    
    % Center coordinate to normalize
    cx = floor(size(x_mesh,1)/2); cy = cx;
    
    % Difference in phase layer 1 and 2: Simulated data
        %delta_ang_sim_l12 = wrapToPi(angle(f_mesh(:,:,1))-angle(f_mesh(:,:,2))- angle(f_mesh(cx,cx+1,1))+angle(f_mesh(cx,cx+1,2))); 
    % Difference in magnitude layer 1 and 2: Simulated data    
        %delta_mag_sim_l12 = abs(f_mesh(:,:,1))-abs(f_mesh(:,:,2));
    % Difference in phase layer 1 and 2: Estimated data
        %delta_ang_est_l12 = wrapToPi(angle(temp(:,:,1))-angle(temp(:,:,2))-angle(temp(cx,cx+1,1))+angle(temp(cx,cx+1,2)));
    % Difference in magnitud layer 1 and 2: Estimated data
        %delta_mag_est_l12 = abs(temp(:,:,2))-abs(temp(:,:,1)); 
    % Layer 1 Simulated/Measured and Layer 1 Estimated
        %delta_mag_sim_est_l1 = abs(f_mesh(:,:,1))-abs(temp(:,:,1)); 
    % Layer 2 Simulated/Measured and Layer 2 Estimated
        %delta_mag_sim_est_l2 = abs(f_mesh(:,:,2))-abs(temp(:,:,2));    
    
    
   % figure('Name','Difference in Magnitude and Phase between layers for Sim. and Est. Values ','units','normalized','outerposition',[0 0 1 1])
   % subplot(2,3,1); cla; surf(x_mesh,y_mesh, delta_mag_sim_l12);  title(['Simulation: | \Delta(L1, L2) |  iter:' num2str(i)]);   colorbar ; shading interp;  view(0,0);
   % subplot(2,3,4); cla; surf(x_mesh,y_mesh, delta_ang_sim_l12);  title(['Simulation: \Delta (\angle L1, \angle L2)  iter:' num2str(i)]);   colorbar ; shading interp;  view(0,90);
   % subplot(2,3,2); cla; surf(x_mesh,y_mesh, delta_mag_est_l12);  title(['Estimation: | \Delta(L1, L2) |  iter: ' num2str(i)]);   colorbar ; shading interp;  view(0,0);
   % subplot(2,3,5); cla; surf(x_mesh,y_mesh, delta_ang_est_l12);  title(['Estimation: \Delta(\angle L1, \angle L2)  iter: ' num2str(i)]);  colorbar ; shading interp;  view(0,90);   
   % subplot(2,3,3); cla; surf(x_mesh,y_mesh, delta_mag_sim_l12-delta_mag_est_l12);  title(['Sim. | \Delta(L1, L2) | - Est. | \Delta(L1, L2) |' num2str(i)]);   colorbar ; shading interp;  view(0,90);
   % subplot(2,3,6); cla; surf(x_mesh,y_mesh, wrapToPi(delta_ang_sim_l12-delta_ang_est_l12));  title(['Sim. \Delta (\angle L1, \angle L2) - Est. \Delta (\angle L1, \angle L2) iter:' num2str(i)]);  colorbar ; shading interp;  view(0,90); 
    
   % figure('Name','Best Case Scenario: error in magnitude ','units','normalized','outerposition',[0 0 1 1])
   % plotIterationsIFT(i*1.0, x_mesh, y_mesh, f_mesh, best_case)
   % subplot(2,4,3);   title(['|Layer 1| -- Best Est. Values iter: ' num2str(iter_best_case(1))]);   colorbar ; shading interp;  view(0,90);
   % subplot(2,4,7);   title(['\angle Layer 1 -- Best Est. Values iter:' num2str(iter_best_case(1))]);   colorbar ; shading interp;  view(0,90);
   % subplot(2,4,4);  title(['|Layer 2| -- Best Est. Values iter: ' num2str(iter_best_case(2))]);   colorbar ; shading interp;  view(0,90);
   % subplot(2,4,8);  title(['\angle Layer 2 -- Best Est. Values iter:' num2str(iter_best_case(2))]);  colorbar ; shading interp;  view(0,90);   

    
   % figure('Name','Best Case Scenario: error in phase ','units','normalized','outerposition',[0 0 1 1])
   % plotIterationsIFT(i*1.0, x_mesh, y_mesh, f_mesh, best_case_p)
   % subplot(2,4,3);   title(['|Layer 1| -- Best Est. Values iter: ' num2str(iter_best_case_p(1))]);   colorbar ; shading interp;  view(0,90);
   % subplot(2,4,7);   title(['\angle Layer 1 -- Best Est. Values iter:' num2str(iter_best_case_p(1))]);   colorbar ; shading interp;  view(0,90);
   % subplot(2,4,4);  title(['|Layer 2| -- Best Est. Values iter: ' num2str(iter_best_case_p(2))]);   colorbar ; shading interp;  view(0,90);
   % subplot(2,4,8);  title(['\angle Layer 2 -- Best Est. Values iter:' num2str(iter_best_case_p(2))]);  colorbar ; shading interp;  view(0,90);   

end 