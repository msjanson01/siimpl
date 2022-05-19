% function plotQTz(z, Q, T);

function plotQTz(z, f1, f2);



% disp(['I(Q + T)dz = ',num2str(sum(f1+f2)*dz),' (keV/cm2)']);

% ------------------------------
f=mjfig;
a=linlogax;
 
plotstairs(z, f1 +eps,'b');
plotstairs(z, f2 +eps,'r');

legend('Q - electronic exitation','T - phonons');
hold on

% adjust min-y on plot (due to log scale and 0+eps!!)
% find smallest fx not equal to 0

win=axis;
slask=[f1 f2];
i=find(slask==0);
slask(i)=max(slask);
axis([win(1:2) 10^floor(log10(min(slask))) 10^ceil(log10(max([f1, f2])))]);
% --------------------------------
title('siimpl - energy depostion')
xlabel('Depth (Å)');
ylabel('Energy loss (keV/Å/cm^{2})');
