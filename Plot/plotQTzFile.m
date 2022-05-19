% plotQTzFile.m


[z, Q, T]=import_mj('QTz.txt');

% disp(['I(Q + T)dz = ',num2str(sum(f1+f2)*dz),' (keV/cm2)']);

% ------------------------------
plotQTz(z, Q, T);