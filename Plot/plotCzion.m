% function plotCzion(z, Cz);

function plotCzion(z, f1);

% ------------------------------
f=mjfig;
a=linlogax;
% -----------------

plotstairs(z, f1 +eps,'b');

mom=moments(z, f1);
legend(['<Rp> =',num2str(mom(1),'%0.0f'),' (Å),  <dRp> =',...
      			  num2str(mom(2),'%0.0f'),' (Å)']);

hold on

% adjust min-y on plot (due to log scale and 0+eps!!)
% find smallest fx not equal to 0

win=axis;
slask=f1;
i=find(slask==0);
slask(i)=max(slask);
axis([win(1:2) 10^floor(log10(min(slask*0.8))) 10^ceil(log10(max([f1*1.2])))]);
% --------------------------------
title('siimpl - ion distribution')
xlabel('Depth (Å)');
ylabel('Czion (cm^{-3})');
