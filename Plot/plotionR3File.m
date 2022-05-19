% plotionR3File.m

A=import_mj('ionR3.txt');

R0 = A(:,1:3);
Rstop = A(:,4:6);


%----------------
f=fig3D;
plotatomsR3(Rstop, 0);

impR = implantarea;

plot3(impR(:,1), impR(:,2), impR(:,3),'-k');
