
function [Px,Py,Pz]=calculatePower(Ex, Ey, Ez, Hx, Hy, Hz)

S = cross([Ex, Ey, Ez],conj([Hx, Hy, Hz])); % Poynting Vector
TAS = real(S);% Time average poynting vector
Px = TAS(:,1);
Py = TAS(:,2);
Pz = TAS(:,3);

end 
%chismo = abs(TaS(:,1)-Px)./abs(Px+0.000000001);
%plot(chismo); 
%legend('calculated', 'sim')


