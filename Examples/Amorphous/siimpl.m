tic
! d:\siimpl\siimpl 115keV_Al_in_aSiC.m
time = toc;

disp(['The simulation completed in ',num2str(time,1),' s']);

plotCzionFile;
load SIMS.mat
plot(xData*1e4, yData,'r');

plotCzIVFile;

