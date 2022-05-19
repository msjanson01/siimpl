% function Se_srim2siimpl(Z1, M1 (amu), Z2, M2 (amu), N, Ef (keV));
% 
% Example: He in SiC
% Z1 = 2; M1 = 4.0; Z2 = [6 14]; M2 = [12 28]; N = [1 1]; Ef = 25 * M1;
% [A1, p, p2, A3, A4, A5] = Se_srim2siimpl(Z1, M1, Z2, M2, N, Ef);
% gives output
%
% A1 =  3.1186   eV/(1015/cm)
% p =   0.5096
% p2 =  0.0405
% A3 =  37.8877 eV/(1015/cm)
% A4 =  2.4099  
% A5 =  4.5551

function [A1, p, p2, A3, A4, A5] = Se_srim2siimpl(Z1, M1, Z2, M2, N, Ef);


% Import SRIM stopping
E_v0 = 25; % keV/amu
Emin = E_v0 * M1  /1000; % = 0.01 of bohr velocity 
Emax = E_v0 * M1 * Z1^(2/3) * 1000; % 1000 * Thomas-Fermi velocity 

E = logspace(log10(Emin), log10(Emax), 100);

[E, Se_srim, Sn_srim] = importSRIMstopping(Z1, M1, Z2, M2, N, E);

if E == 0
    A1 = 0;
    p = 0;
    p2 = 0;
    A2 = 0;
    A4 = 0;
    A5 = 0;
    disp('ERROR!');
    return;
end

% fit Se_siimpl stopping;

% Start guess
[A1, p, p2, A3, A4, A5] = Se_siimpl_default(Z1, M1, Z2, N);

[A1, p, p2, A3, A4, A5] = fit_Se_siimpl(E , Se_srim, Ef, A1, p, p2, A3, A4, A5);

