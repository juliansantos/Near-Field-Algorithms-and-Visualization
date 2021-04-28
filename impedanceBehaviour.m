% This function plot the behaviour of the impedance at different layers,
% this is with the objective to see until what point the assumption of
% 377ohm is valid for the caracterization of the antenna. 


% ToDo: plot the surface of impedance in one layer

% ToDo: plot the mean impedance with std intervals with respect to z

function impedanceBehaviour(X,Y,Z,Ex, Ey, Ez, Hx, Hy, Hz,layers)
    zref = 41.2 ; % Position in z direction of the aperture
    eta = sqrt(abs(Ex).^2+abs(Ey).^2+abs(Ez).^2)./sqrt(abs(Hx).^2+abs(Hy).^2+abs(Hz).^2);
    eta_layers = reshape(eta,[length(unique(X)),length(unique(Y)),length(unique(Z))]);
    u_eta= squeeze(mean(mean(eta_layers))); 
    s_eta = squeeze(std(std(eta_layers))); 
    
    % Uncertainity plot
    subplot(2,2,[1,2]);
    errorbar(unique(Z)-zref,u_eta,s_eta,'-s','MarkerSize',7,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')
    grid on;
    xlim([-1 53]); ylim([347 415]);
    xlabel('Distance from the aperture [mm]')
    ylabel('Impedance [\Omega]')
    legend('\mu, \sigma')
    title('Impedance behaviour in z direction')
    
    %plot(eta(112361:end)) % plot for all x and y
    [val,idx]=min(abs(layers(1) + zref - Z)); z(1)= Z(idx);
    [val,idx]=min(abs(layers(2) + zref - Z)); z(2)= Z(idx);
    
    xtemp = X(Z==z(1));
    x_mesh = reshape(xtemp, [length(unique(X)),length(unique(Y))]);
    ytemp = Y(Z==z(1));
    y_mesh = reshape(ytemp, [length(unique(Y)), length(unique(X))]);
    e_temp = eta(Z==z(1));
    e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
    subplot(2,2,3);
    surf(x_mesh,y_mesh,e_mesh)
    colorbar
    xlabel('x [mm]')
    ylabel('y [mm]')
    title(['Impedance at ' num2str(z(1)-zref) 'mm from the aperture']) 
    view(0,90)
    
    subplot(2,2,4);
    e_temp = eta(Z==z(2));
    e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
    surf(x_mesh,y_mesh,e_mesh)
    colorbar
    xlabel('x [mm]')
    ylabel('y [mm]')
    title(['Impedance at ' num2str(z(2)-zref) 'mm from the aperture']) 
    view(0,90)
    % make xmesh and y mesh and get plot the surface of the behaviour of
    % the impedance
    
end 