E = lindist(1040, 1150, 1000);

nu = 180;
alpha = nu;

dE_det = 1;

yield = siimpl_RBS(E, dE_det, nu, alpha);

mjaxes;
plot(E, yield);
xlabel('Energy (keV)');
ylabel('Yield (a.u.)');
axis([min(E) max(E) 0 0.5])