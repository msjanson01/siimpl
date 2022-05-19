% plotionR3File.m

A=import_mj('ionR3.txt');

R0 = A(:,1:3);
Rstop = A(:,4:6);

Rrel=Rstop-R0;

%----------------
f=fig3D;
plotatomsR3(Rrel, 0);

R2 = max(sum(Rrel.^2,2));

impDir = implantdir;


Rdir = sqrt(R2)*impDir;

plot3([0 Rdir(1)], [0 Rdir(2)], [0 Rdir(3)], 'k'); 


