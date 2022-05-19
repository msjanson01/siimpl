E = lindist(1100, 1150, 100);

nu = 170;
alpha = nu;

dE_det = 1;

yield = siimpl_RBS(E, dE_det, nu, alpha);

mjaxes;
plot(E, yield);
