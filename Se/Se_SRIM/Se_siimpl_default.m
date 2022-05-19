% function [A1, p, p2, A3, A4, A5] = Se_siimpl_default(Z1, M1, Z2, N);

function [A1, p, p2, A3, A4, A5] = Se_siimpl_default(Z1, M1, Z2, N);

% Se_lo
A1 =  3.85454e-2 / sqrt(M1) * Z1.^(7/6) * Z2 ./ (Z1.^(2/3) + Z2.^(2/3)).^1.5 * 1e2; % eV/(1015cm-2) cm^2
A1 = sum(A1 .* N)/sum(N);

p = 0.5;
p2 = 0.0;

% Se_hi
Z2eff = sum(Z2 .* N)/sum(N);

if (Z2eff >= 13)
   I0 = (9.76 + 58.5 * Z2eff^-1.19); 
else
   I0 = (12 + 7/Z2eff); 
end  

I = I0 * Z2eff;
 
if (Z1 < 3)
   C = 100 * Z1/Z2eff;
else
   C = 5;
end   

A3 = 0.238 * Z1^2 * Z2eff * M1; %   keV eV/(1015cm-2)
A4 = 4.55e-4 * C * I * M1; %	    keV				(Xb)
A5 = 2.20e3 / (I * M1); %        keV-1				(Xc)

