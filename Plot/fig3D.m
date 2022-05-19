% function fig3D = fig3D;

function fig3D = fig3D;


% init 3D fig and axes ------------------------------
fig3D=figure;
set(fig3D,'color','w');

ax3D=axes;
plot3(0,0,0);
set(ax3D,'color','w');
cla;
box on;
axis equal;
axis on;
hold on;
set(ax3D,'proj','ortho');
xlabel('X');
ylabel('Y');
zlabel('Z');



set(fig3D,'PaperUnits','centimeter','PaperOrientation','landscape');
set(fig3D,'PaperType','a4');
set(fig3D,'PaperPosition',[3 2 24 16.5]);
