% Antenna Diagnostics - Back Projection using Planar Near-Field
% Measurements
% AUTHOR: B*log10(Engineering)
% WEBSITE: http://blog.p86.dk/
% DATE: 25 Jan 2012
%
% DESCRIPTION:
% This script is an example implemenation of back projecting using planar
% near-field measurements.

% Clear Workspace
clear
clc
close all

% Formatting Variables
format_fontsize_title = 18;
format_fontsize_axis = 16;

% Load Measurement Data
% Standard Gain Horn, full span 
% x = 0:0.025:1.5 [m]
% y = 0:0.025:0.8 [m]
% Contain measurements for 5.2, 5.25 and 5.3 GHz
load('SGH_FullPlanarScan.mat');

f = 5.3e9; % Frequency of Interest [Hz]
z_0 = 0.088; % Distance between the AUT and probe apertures [m]
delta_x = 0.022; % planar x scan step size [m]
delta_y = 0.022; % planar y scan step size [m]

% Perform Back Projection to the Aperture of the AUT
[E_x E_y E_z, x, y] = BackProjection_PlanarNearField(comp_x_measurement(:,:,3),comp_y_measurement(:,:,3),z_0,f,delta_x,delta_y);

% Rotate Axis such that one views the AUT from the front
E_x = abs(fliplr(E_x'));
E_y = abs(fliplr(E_y'));
E_z = abs(fliplr(E_z'));

% Normalise Fields
E_x_norm = E_x./max(max(E_x));
E_y_norm = E_y./max(max(E_y));
E_z_norm = E_z./max(max(E_z));

% Convert to dB scale
E_x_mag = 20*log10(abs(E_x_norm));
E_y_mag = 20*log10(abs(E_y_norm));
E_z_mag = 20*log10(abs(E_z_norm));

% Find Power Intensity
U = abs(E_x).^2 + abs(E_y).^2 + abs(E_z).^2;
U_dB = 10*log10(U); % Convert to dB scale

% Plot Data
figure('Name','U - 3.95-5.85 GHz Standard Gain Horn Back Projection')
set(gca,'fontsize',format_fontsize_axis)
surf(x*1000,y*1000,U_dB);
title({'3.95-5.85 GHz Standard Gain Horn Back Projection','Frequency: 5.3 GHz , z_0=88 mm'},'fontsize',format_fontsize_title)
xlabel('x (mm)','fontsize',format_fontsize_axis);
ylabel('y (mm)','fontsize',format_fontsize_axis);
zlabel('U_{norm} (dB)','fontsize',format_fontsize_axis);
set(gca,'XLim',[min(x)*1000 max(x)*1000]);
set(gca,'YLim',[min(y)*1000 max(y)*1000]);
view(-37.5,30);
shading flat;
colorbar;
daspect([1 1 1]);

figure('Name','Ex - 3.95-5.85 GHz Standard Gain Horn Back Projection')
set(gca,'fontsize',format_fontsize_axis)
surf(x*1000,y*1000,E_x_mag);
title({'3.95-5.85 GHz Standard Gain Horn Back Projection','Frequency: 5.3 GHz , z_0=88 mm'},'fontsize',format_fontsize_title)
xlabel('x (mm)','fontsize',format_fontsize_axis);
ylabel('y (mm)','fontsize',format_fontsize_axis);
zlabel('U_{norm} (dB)','fontsize',format_fontsize_axis);
set(gca,'XLim',[min(x)*1000 max(x)*1000]);
set(gca,'YLim',[min(y)*1000 max(y)*1000]);
view(-37.5,30);
shading flat;
colorbar;
daspect([1 1 1]);

figure('Name','Ey - 3.95-5.85 GHz Standard Gain Horn Back Projection')
set(gca,'fontsize',format_fontsize_axis)
surf(x*1000,y*1000,E_y_mag);
title({'3.95-5.85 GHz Standard Gain Horn Back Projection','Frequency: 5.3 GHz , z_0=88 mm'},'fontsize',format_fontsize_title)
xlabel('x (mm)','fontsize',format_fontsize_axis);
ylabel('y (mm)','fontsize',format_fontsize_axis);
zlabel('U_{norm} (dB)','fontsize',format_fontsize_axis);
set(gca,'XLim',[min(x)*1000 max(x)*1000]);
set(gca,'YLim',[min(y)*1000 max(y)*1000]);
view(-37.5,30);
shading flat;
colorbar;
daspect([1 1 1]);

figure('Name','Ez - 3.95-5.85 GHz Standard Gain Horn Back Projection')
set(gca,'fontsize',format_fontsize_axis)
surf(x*1000,y*1000,E_z_mag);
title({'3.95-5.85 GHz Standard Gain Horn Back Projection','Frequency: 5.3 GHz , z_0=88 mm'},'fontsize',format_fontsize_title)
xlabel('x (mm)','fontsize',format_fontsize_axis);
ylabel('y (mm)','fontsize',format_fontsize_axis);
zlabel('U_{norm} (dB)','fontsize',format_fontsize_axis);
set(gca,'XLim',[min(x)*1000 max(x)*1000]);
set(gca,'YLim',[min(y)*1000 max(y)*1000]);
view(-37.5,30);
shading flat;
colorbar;
daspect([1 1 1]);


%% Test with CST 10GHz antenna

clc; clear all; close all force; 

load DataCST.mat

f = 10e9; % Frequency of Interest [Hz]
Field = Ex_40mm;
Desired_Field = Ex_0mm; 
z_0 = 0.040; % Distance between the AUT and probe apertures [m]
delta_x = 0.001; % planar x scan step size [m]
delta_y = 0.001; % planar y scan step size [m]

[E_x, x, y] = BackProjection_PlanarNearField_v2(Field,z_0,f,delta_x,delta_y);

      figure; 
      subplot(1,3,1)
        surf(x_mesh, y_mesh, abs(Field)); shading interp; colorbar; view(90,90)
        title(['Input - Max value = ' num2str(max(max(abs(Field))))])
      subplot(1,3,2)
        surf(x_mesh, y_mesh, abs(E_x)); shading interp; colorbar; view(90,90)
        title(['Output - Est. Mag' num2str(max(max(abs(E_x))))])
      subplot(1,3,3)
        surf(x_mesh, y_mesh, abs(Desired_Field)); shading interp; colorbar; view(90,90)
        title(['Desired. Mag. - Max value = ' num2str(max(max(abs(Desired_Field))))])
        
      figure; 
      subplot(1,3,1)
        surf(x_mesh, y_mesh, abs(Field)); shading interp; colorbar; view(180,0)
        title(['Input - Max value = ' num2str(max(max(abs(Field))))])
      subplot(1,3,2)
        surf(x_mesh, y_mesh, abs(E_x)); shading interp; colorbar; view(180,0)
        title(['Output - Est. Mag' num2str(max(max(abs(E_x))))])
      subplot(1,3,3)
        surf(x_mesh, y_mesh, abs(Desired_Field)); shading interp; colorbar; view(180,0)
        title(['Desired. Mag. - Max value = ' num2str(max(max(abs(Desired_Field))))]) 
        
        
    figure; 
      subplot(1,3,1)
        surf(x_mesh, y_mesh, angle(Field)); shading interp; colorbar; view(90,90)
        title(['Input '])
      subplot(1,3,2)
        surf(x_mesh, y_mesh, angle(E_x)); shading interp; colorbar; view(90,90)
        title(['Output '])
      subplot(1,3,3)
        surf(x_mesh, y_mesh, angle(Desired_Field)); shading interp; colorbar; view(90,90)
        title(['Desired Phase'])