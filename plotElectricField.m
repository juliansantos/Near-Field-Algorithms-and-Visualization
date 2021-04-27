function plotElectricField(X, Y, Z, Ex, Ey,Ez)
    figure;
    subplot(2,4,[1,2,5,6])
    quiver3(X, Y, Z, abs(Ex), abs(Ey), abs(Ez),'LineWidth',4,'color','b');
    title('Electrical Field - Magnitud - 60GHz Horn Antenna')
    axis equal
    hold on; plotstructure3D()
    % 43.2 46.2 51.2
    z = [51.1 20.1]; % vector with the alpha cuts
    
    % TODO: check length of z, to verify is two. and existence of those
    % values. input argument, pos in z of aperture. 
    % ===convert the value z to the real distance from the aperture and not
    % the distance from the reference position. 
    
    subplot(2,4,4)
    xtemp = X(find(Z==z(1)));
    x_mesh = reshape(xtemp, [length(unique(X)), length(xtemp)/length(unique(X))]);
    ytemp = Y(find(Z==z(1)));
    y_mesh = reshape(ytemp, [length(unique(Y)), length(ytemp)/length(unique(Y))]);
    e_temp = Ex(find(Z==z(1)));
    e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
    pcolor(x_mesh,y_mesh,((angle(e_mesh))))
    shading interp; 
    colorbar; colormap(jet);
    xlabel('x [mm]')
    ylabel('y [mm]')
    title(['Phase of Ex at ' num2str(z(1)) 'mm distance from the aperture']) 
    
    subplot(2,4,3)
    xtemp = X(find(Z==z(1)));
    x_mesh = reshape(xtemp, [length(unique(X)), length(xtemp)/length(unique(X))]);
    ytemp = Y(find(Z==z(1)));
    y_mesh = reshape(ytemp, [length(unique(Y)), length(ytemp)/length(unique(Y))]);
    e_temp = Ex(find(Z==z(1)));
    e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
    pcolor(x_mesh,y_mesh,(mag2db(abs(e_mesh))))
    shading interp; 
    colorbar; colormap(jet);
    xlabel('x [mm]')
    ylabel('y [mm]')
    title(['Magnitude of Ex at ' num2str(z(1)) 'mm distance from the aperture']) 

    subplot(2,4,8)
    xtemp = X(find(Z==z(2)));
    x_mesh = reshape(xtemp, [length(unique(X)), length(xtemp)/length(unique(X))]);
    ytemp = Y(find(Z==z(2)));
    y_mesh = reshape(ytemp, [length(unique(Y)), length(ytemp)/length(unique(Y))]);
    e_temp = Ex(find(Z==z(2)));
    e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
    pcolor(x_mesh,y_mesh,((angle(e_mesh))))
    shading interp; 
    colorbar; colormap(jet);
    xlabel('x [mm]')
    ylabel('y [mm]')
    title(['Phase of Ex at ' num2str(z(2)) 'mm distance from the aperture']) 
    
    subplot(2,4,7)
    xtemp = X(find(Z==z(2)));
    x_mesh = reshape(xtemp, [length(unique(X)), length(xtemp)/length(unique(X))]);
    ytemp = Y(find(Z==z(2)));
    y_mesh = reshape(ytemp, [length(unique(Y)), length(ytemp)/length(unique(Y))]);
    e_temp = Ex(find(Z==z(2)));
    e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
    pcolor(x_mesh,y_mesh,(mag2db(abs(e_mesh))))
    shading interp; 
    colorbar; colormap(jet);
    xlabel('x [mm]')
    ylabel('y [mm]')
    title(['Magnitude of Ex at ' num2str(z(2)) 'mm distance from the aperture'])
end 