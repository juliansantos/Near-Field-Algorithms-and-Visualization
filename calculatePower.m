
clc;

S = cross([Ex, Ey, Ez],conj([Hx, Hy, Hz])); % Poynting Vector
W  = 0.5*abs([Ex, Ey, Ez])+0.5*abs([Hx, Hy, Hz]);


dif = [real(S)./Px];

plot(dif(:,1))