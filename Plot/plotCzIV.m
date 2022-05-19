% function plotIVz(z, Q, T);

function plotIVz(z, f1, f2);

[m,n] = size(f1);

Icol={'b','r','k','g','m','c','y'};
Vcol={':b',':r',':k',':g',':m',':c',':y'};

% ------------------------------
FIG=mjfig;
AX=linlogax;

f=[];
for i=1:m
	plotstairs(z, f1(i,:) +eps, Icol{i});
	plotstairs(z, f2(i,:) +eps, Vcol{i});
	f=[f,f1(i,:),f2(i,:)];    
end

legendcell={};
for i=1:m
   legendcell={legendcell{:}, ['I ',num2str(i)], ['V ',num2str(i)]};
end


legend(legendcell);
hold on
title('siimpl - point defects')
xlabel('Depth (Å)');
ylabel('Point defects (cm^{-3})');

if max(f) == 0, return, end

% adjust min-y on plot (due to log scale and 0+eps!!)
% find smallest fx not equal to 0
% adjust min-y on plot (due to log scale and 0+eps!!)
% find smallest fx not equal to 0
win=axis;
slask=f;
i=find(slask==0);
slask(i)=max(slask);

axis([win(1:2) 10^floor(log10(min(slask))) 10^ceil(log10(max([f])))]);
% --------------------------------

