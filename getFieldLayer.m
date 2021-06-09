% Function to get the layer at a distance x in var. layers from the aperture. 
%       This function assumes that all layers have the same number of
%       discrete elements over XcrossY

function [x_mesh, y_mesh, f_mesh]= getFieldLayer(X, Y, Z, F, layers)
    z_values = unique(Z);
    flag = 0; 
    xref = 0; yref = 0; zref = 41.2 ; % Reference of the aperture of the antenna
    layers = layers + zref;
    x_mesh=[]; y_mesh=[]; f_mesh=[];
    for j=1:length(layers)
        flag = flag + sum(layers(j)==z_values);
    end
    
    if flag ~= length(layers)
        fprintf(2,'There are one or more values of the selected layers that are not in the introduced dataset. ')       
    else    
        Lux = length(unique(X)); % Length Unique Elements at X
        Luy = length(unique(Y)); % Length Unique Elements at Y
        xtemp = X(Z==layers(1));
        x_mesh = reshape(xtemp, [Lux, length(xtemp)/Lux]);
        %ytemp = Y(Z==layers(1));
        %y_mesh = reshape(ytemp, [Luy, length(ytemp)/Luy]);
        y_mesh = x_mesh';
        dimF = [size(x_mesh,1),size(y_mesh,2)] ;  
        
        f_mesh = zeros([dimF,length(layers)]);
        for i=1:length(layers)
            f_temp = F(Z==layers(i));
            f_mesh(:,:,i) = reshape(f_temp, dimF);
        end 
    end
end 