% plotIVR3File.m

A=import_mj('VR3.txt');
Vid = A(:,1);
RV = A(:,2:4);

A=import_mj('IR3.txt');
Iid = A(:,1);
RI = A(:,2:4);

%----------------
f=fig3D;
box on;
title('siimpl - self interstitials and vacancies')

plotatomsR3(RI, Iid);
plotatomsR3(RV, -Vid);

impR = implantarea;

plot3(impR(:,1), impR(:,2), impR(:,3),'-k');

