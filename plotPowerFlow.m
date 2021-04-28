function plotPowerFlow(X, Y, Z, Ex, Ey, Ez, layers, size_layers)
    % Input parameters: 
    %                   X
    
    xref = 0; yref = 0; zref = 41.2 ; %Describe the position with respect to [0,0,0] at 
        % which is located the center of the aperture of the horn antenna. 
    z = [0,0];

    figure;
    
    switch nargin 
        case 6 
            subplot(2,4,1:8)
        case 8
            [val,idx]=min(abs(layers(1) + zref - Z)); z(1)= Z(idx);
            [val,idx]=min(abs(layers(2) + zref - Z)); z(2)= Z(idx);
            subplot(2,5,[1,2,6,7])
        otherwise 
    end 

    quiver3(X, Y, Z, Ex, Ey, Ez,'LineWidth',4,'color',[0.5 0 0.5]);
    title('PowerFlow - Magnitude - 60GHz Horn Antenna')
    axis equal
    hold on; 
    
    switch nargin 
        case 6 
            plotstructure3D()
        case 8
            plotstructure3D(z, size_layers, [0 0 1])
           

            subplot(2,5,3)
            xtemp = X(Z==z(1));
            x_mesh = reshape(xtemp, [length(unique(X)), length(unique(Y))]);
            ytemp = Y(Z==z(1));
            y_mesh = reshape(ytemp, [length(unique(Y)), length(unique(X))]);
            e_temp = Ex(Z==z(1));
            e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
            pcolor(x_mesh,y_mesh,(mag2db(abs(e_mesh))))
            shading interp; 
            colorbar; colormap(jet);
            xlabel('x [mm]')
            ylabel('y [mm]')
            title(['Magnitude of Px at ' num2str(z(1)-zref) 'mm from the AUT']) 

            subplot(2,5,4)
            e_temp = Ey(Z==z(1));
            e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
            pcolor(x_mesh,y_mesh,(mag2db(abs(e_mesh))))
            shading interp; 
            colorbar; colormap(jet);
            xlabel('x [mm]')
            ylabel('y [mm]')
            title(['Magnitude of Py at ' num2str(z(1)-zref) 'mm from the AUT']) 

            subplot(2,5,5)
            e_temp = Ez(Z==z(1));
            e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
            pcolor(x_mesh,y_mesh,(mag2db(abs(e_mesh))))
            shading interp; 
            colorbar; colormap(jet);
            xlabel('x [mm]')
            ylabel('y [mm]')
            title(['Magnitude of Pz at ' num2str(z(1)-zref) 'mm from the AUT']) 
            
            subplot(2,5,8)
            e_temp = Ex(Z==z(2));
            e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
            pcolor(x_mesh,y_mesh,(mag2db(abs(e_mesh))))
            shading interp; 
            colorbar; colormap(jet);
            xlabel('x [mm]')
            ylabel('y [mm]')
            title(['Magnitude of Px at ' num2str(z(2)-zref) 'mm from the AUT'])
            
            subplot(2,5,9)
            e_temp = Ey(Z==z(2));
            e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
            pcolor(x_mesh,y_mesh,(mag2db(abs(e_mesh))))
            shading interp; 
            colorbar; colormap(jet);
            xlabel('x [mm]')
            ylabel('y [mm]')
            title(['Magnitude of Py at ' num2str(z(2)-zref) 'mm from the AUT'])
            
            subplot(2,5,10)
            e_temp = Ez(Z==z(2));
            e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
            pcolor(x_mesh,y_mesh,(mag2db(abs(e_mesh))))
            shading interp; 
            colorbar; colormap(jet);
            xlabel('x [mm]')
            ylabel('y [mm]')
            title(['Magnitude of Pz at ' num2str(z(2)-zref) 'mm from the AUT'])
        otherwise 
    end    
end 

