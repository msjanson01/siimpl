% function yield = siimpl_RBS(E (keV), dE_det (keV), nu (deg), alpha (deg));
function yield = siimpl_RBS(E, dE_det, nu, alpha);

nu = nu*pi/180;
alpha = alpha*pi/180;


[Z2 M2 Natom] = target_atoms;
Natom = Natom * 1e-24; % Å^-3

[Z1, M1, E0] = implant_ion;

impdir = implantdir;
snorm = siimpl_surfnorm;

theta = asin( norm(cross(impdir, -snorm)));

[Efermi, A1, A2, p2, A3, A4, A5] = Se_parameters;

A=import_mj('NEP.txt');

if isempty(A),disp('Can not open NEP.txt file!'); return; end

[m,n] = size(A);
z = A(:,1)';
NEP = A(:,2:n)';

yield = zeros(size(E));

dz = z(2) - z(1);


Se_fun = inline('Se_siimpl(A, E_pp,  Efermi)*1e-2','E_pp','A','Efermi');


for i = 1:length(Z2)
    
    % E_p = energy of beam in
    E_p = E0;
    
    for j = 1:length(z)
        
        % --- ion out
        % E_pp = energy of beam out
        
        % Rutherford scattering
        E_pp = (1 - 4*M1*M2(i)/(M1 + M2(i))^2*sin(nu/2)^2 ) * E_p;
        

        % Constant Se for out trajectory, updated later
        Se = Se_siimpl([A1 A2 p2 A3 A4 A5], E_pp,  Efermi); % eV/(1e15 atoms/cm-2)
        Se = Se *1e-2;  %  eV/(10^15cm^-2) -> keV*Å^2
        Eout = E_pp - Natom * Se * z(i)/(cos(pi - alpha));

        % Calculate yield
        dE = dE_det; % Add dQ later!!!
        yield = yield + NEP(i,j)/sqrt(2*pi)/dE * exp(-((E - Eout)/dE).^2/2);
        
        % --- Update energy of beam
        Se = Se_siimpl([A1 A2 p2 A3 A4 A5], E_p,  Efermi); % eV/(1e15 atoms/cm-2)
        Se = Se *1e-2;  %  eV/(10^15cm^-2) -> keV*Å^2
        E_p = E_p - Natom * Se * dz/cos(theta);
        
    end
  
end