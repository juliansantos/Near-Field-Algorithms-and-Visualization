% This code calculate the Far Filed pattern 

function calculateFarField 
 n = 0; % Impedance in free Space
 lambda = 5e-3;
 k = 0; % Wave number

 theta = 0; 
 Er = n*cos(theta)(1+1/(1j*k*r))*exp(-1j*k*r)/(2*pi*r^2);
 Ephi = 1j*n*cos(theta)*(1+1/(1j*k*r) -1/(1j*k^2*r^2) )*exp(-1j*k*r)/(4*pi*r);
 Eth = 0;
 
end
 
 