% Author: Julián Santos
% Description: This code plot an easy representation of the 3D horn antenna
% 60GHz
%% Plot 3D structure and coordinate system
function plotstructure3D
vert = [-5.1 -6.565 41.2;5.1 -6.565 41.2;5.1 6.565 41.2;-5.1 6.565 41.2];
patch('Vertices',vert,'Faces',[1 2 3 4],'FaceVertexCData',hsv(2),'FaceColor',[1 0.84 0])

vert = [-3.6 -5.065 41.2;3.6 -5.065 41.2;3.6 5.065 41.2;-3.6 5.065 41.2];
patch('Vertices',vert,'Faces',[1 2 3 4],'FaceVertexCData',hsv(2),'FaceColor',[1 0 0])

vert = [-2.39 -3.2 13;2.39 -3.2 13;2.39 3.2 13;-2.39 3.2 13];
patch('Vertices',vert,'Faces',[1 2 3 4],'FaceVertexCData',hsv(2),'FaceColor',[1 0.84 0])

vert = [-0.89 -1.7 13;0.89 -1.7 13;0.89 1.7 13;-0.89 1.7 13];
patch('Vertices',vert,'Faces',[1 2 3 4],'FaceVertexCData',hsv(2),'FaceColor',[1 0 0])

line([-2.39 -2.39],[3.2 3.2],[13 26],'Color','black','LineStyle','-')
line([2.39 2.39],[3.2 3.2],[13 26],'Color','black','LineStyle','-')
line([-2.39 -2.39],[-3.2 -3.2],[13 26],'Color','black','LineStyle','-')
line([2.39 2.39],[-3.2 -3.2],[13 26],'Color','black','LineStyle','-')


line([-2.39 -5.1],[3.2 6.565],[26 41.2],'Color','black','LineStyle','-')
line([2.39 5.1],[3.2 6.565],[26 41.2],'Color','black','LineStyle','-')
line([-2.39 -5.1],[-3.2 -6.565],[26 41.2],'Color','black','LineStyle','-')
line([2.39 5.1],[-3.2 -6.565],[26 41.2],'Color','black','LineStyle','-')

hold on;
% Plot coordinate system
quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[10;0;0],[0;10;0],[0;0;10],'color','black'); 

beta = 1;
L=50; %length of the layer of measurement in milimeters
z = [51.1 20.1];
if beta == 1 
    vert = [-L/2 -L/2 z(1);L/2 -L/2 z(1);L/2 L/2 z(1);-L/2 L/2 z(1)];
    patch('Vertices',vert,'Faces',[1 2 3 4],'FaceVertexCData',hsv(2),'FaceColor',[1 0 0])
    
    vert = [-L/2 -L/2 z(2);L/2 -L/2 z(2);L/2 L/2 z(2);-L/2 L/2 z(2)];
    patch('Vertices',vert,'Faces',[1 2 3 4],'FaceVertexCData',hsv(2),'FaceColor',[1 0 0])
end

alpha(.5)
view(21,14) % default view, [az, el]
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')


    

end