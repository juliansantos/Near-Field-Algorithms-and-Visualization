%[MagE, Theta, Phi] Author: J. A. Santos Bustos

function plotFarField(FFdata1,FFdata2,FFdata3)
% Far field 'Normal way' [MagE, Theta, Phi] 3D
    patternCustom(FFdata1,FFdata2,FFdata3);
    title('Far Field 60GHz Horn Antenna');

% Far field 'Rectangular way' 3D
    figure
    patternCustom(FFdata1,FFdata2,FFdata3,...
    'CoordinateSystem','rectangular');
    title('Far Field 60GHz Horn Antenna -- Rectangular Representation');

% 2D Slice 
    figure;
    title('Far Field 60GHz Horn Antenna -- Slice Representation');
    subplot(1,2,1);
%    patternCustom(FFdata1,FFdata2,FFdata3,'CoordinateSystem','rectangular',...
%    'Slice','phi','SliceValue',90);
    subplot(1,2,2);
%    patternCustom(FFdata1,FFdata2,FFdata3,'CoordinateSystem','polar',...
%    'Slice','phi','SliceValue',90);
end 