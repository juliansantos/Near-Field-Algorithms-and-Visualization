function plotElectricField(X, Y, Z, Ex, Ey, Ez, layers, size_layers)
    % Input parameters: 
    %                   X
    
    xref = 0; yref = 0; zref = 167.5 ; %Describe the position with respect to [0,0,0] at 
        % which is located the center of the aperture of the horn antenna. 
    z = [0,0];

    figure('Name','Electric Field on two plannar regions','units','normalized','outerposition',[0 0 1 1])
    
    switch nargin 
        case 6 
            subplot(2,4,1:8)
        case 8
            [val,idx]=min(abs(layers(1) + zref - Z)); z(1)= Z(idx);
            [val,idx]=min(abs(layers(2) + zref - Z)); z(2)= Z(idx);
            subplot(2,4,[1,2,5,6])
        otherwise 
    end 

    quiver3(X, Y, Z, abs(Ex), abs(Ey), abs(Ez),'LineWidth',4,'color','b');
    title('Electrical Field - Magnitude - 60GHz Horn Antenna')
    axis equal
    hold on; 
    
    switch nargin 
        case 6 
            plotstructure3D()
        case 8
            plotstructure3D(z, size_layers, [1 0 0])
            subplot(2,4,4)
            xtemp = X(Z==z(1));
            x_mesh = reshape(xtemp, [length(unique(X)), length(unique(Y))]);
           % ytemp = x_mesh;
            ytemp = Y(Z==z(1));
            y_mesh = reshape(ytemp, [length(unique(X)), length(unique(Y))]);
           % y_mesh = x_mesh';
            e_temp = Ex(Z==z(1));
            e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
            pcolor(x_mesh,y_mesh,((angle(e_mesh))))
            shading interp; 
            colorbar; colormap(jet);
            xlabel('x [mm]')
            ylabel('y [mm]')
            title(['Phase of Ex at ' num2str(z(1)-zref) 'mm distance from the aperture']) 

            subplot(2,4,3)
            e_temp = Ex(Z==z(1));
            e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
            pcolor(x_mesh,y_mesh,(mag2db(abs(e_mesh))))
            shading interp; 
            colorbar; colormap(jet);
            xlabel('x [mm]')
            ylabel('y [mm]')
            title(['Magnitude of Ex at ' num2str(z(1)-zref) 'mm distance from the aperture']) 

            subplot(2,4,8)
            e_temp = Ex(Z==z(2));
            e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
            pcolor(x_mesh,y_mesh,((angle(e_mesh))))
            shading interp; 
            colorbar; colormap(jet);
            xlabel('x [mm]')
            ylabel('y [mm]')
            title(['Phase of Ex at ' num2str(z(2)-zref) 'mm distance from the aperture']) 

            subplot(2,4,7)
            e_temp = Ex(Z==z(2));
            e_mesh = reshape(e_temp, [size(x_mesh,1),size(y_mesh,2)]);
            pcolor(x_mesh,y_mesh,(mag2db(abs(e_mesh))))
            shading interp; 
            colorbar; colormap(jet);
            xlabel('x [mm]')
            ylabel('y [mm]')
            title(['Magnitude of Ex at ' num2str(z(2)-zref) 'mm distance from the aperture'])
        otherwise 
    end    
end 