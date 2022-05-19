tic
! d:\siimpl\siimpl 60keV_Al_in_SiC.m
time = toc;

disp(['The simulation completed in ',num2str(time,2),' s']);

plotCzionFile;
plotCzIVFile;
plotQTzFile;

