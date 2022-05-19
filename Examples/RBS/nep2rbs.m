% NOTE that both the rtelevant NEP.txt and siimplLog.txt files 
% have to be present in the same directory as this file


% Define Energy spectra
Emin = 1100; % (keV)
Emax = 1150; % (keV)
dE = 1; % (keV), width of one channel

E = linspace(Emin, Emax, (Emax - Emin)/dE);


% Detector geometry

% Scatter angle nu (beam - detector) 
nu = 170;
% Angel between surface normal and detector
alpha = nu;

% Detectoe energy resolution
dE_det = 2; % (keV)

% Calculate RBS yield at detected energy E
yield = siimpl_RBS(E, dE_det, nu, alpha);

figure;
axes;

xlabel('Energy (keV)');
ylabel('Backscattering Yield')

plot(E, yield,'-k');
