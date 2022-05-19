% SRIMstopping

% Define ion
Z1 = 8;
M1 = 16;

% Define target
Z2 = [6 14];
M2 = [12 28];
N =  [1 1];

% Import SRIM stopping
E_v0 = 25; % keV/amu
Emin = E_v0 * M1  /1000; % = 0.01 of bohr velocity 
Emax = E_v0 * M1 * Z1^(2/3) * 1000; % 1000 * Thomas-Fermi velocity 

E = expdist(Emin, Emax, 100);

[E, Se_srim, Sn_srim] = importSRIMstopping(Z1, M1, Z2, M2, N, E);
if E == 0
    return
end


figure;
a = loglogax;
xlabel('Energy (keV)');
ylabel('Stopping [eV/(10^{15}/cm^{2})]');

l = plot(E, Se_srim,'k.');


% fit Se_siimpl stopping;
Ef = 25 * M1; % keV

% Start guess
[A1, p, p2, A3, A4, A5] = Se_siimpl_default(Z1, M1, Z2, N);
%plot(E, Se_siimpl([A1, p, p2, A3, A4, A5], E,  Ef), 'k:');

[A1, p, p2, A3, A4, A5] = fit_Se_siimpl(E , Se_srim, Ef, A1, p, p2, A3, A4, A5)

plot(E, Se_siimpl([A1, p, p2, A3, A4, A5], E,  Ef), 'k-');

% ----------


%h = legend('Se_{SRIM}','Sn_{SRIM}');



