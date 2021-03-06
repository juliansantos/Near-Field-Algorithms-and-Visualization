% The input parameter 'beta' allow to select between different types of
% Structure Plot:
%               1. No input -> Only the 3D structure + coordinate system
%               2. Beta = 1 -> The 3D structure + 2 layers 

function plotstructure3D(z, sl, color)
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

    hold on; alpha(.5);
    % Plot coordinate system
    quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[10;0;0],[0;10;0],[0;0;10],'color','black'); 
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')             
    switch nargin
            case 0 
            case 3
                L11 =sl(1,1); L12=sl(1,2); L21=sl(2,1); L22=sl(2,2); 
                vert = [-L11/2 -L12/2 z(1);L11/2 -L12/2 z(1);L11/2 L12/2 z(1);-L11/2 L12/2 z(1)];
                patch('Vertices',vert,'Faces',[1 2 3 4],'FaceVertexCData',hsv(2),'FaceColor',color)

                vert = [-L21/2 -L22/2 z(2);L21/2 -L22/2 z(2);L21/2 L22/2 z(2);-L21/2 L22/2 z(2)];
                patch('Vertices',vert,'Faces',[1 2 3 4],'FaceVertexCData',hsv(2),'FaceColor',color)
                alpha(.5)
                view(21,14) % default view, [az, el]
            otherwise
                % :v
     end
end