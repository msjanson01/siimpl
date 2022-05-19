% function plotIVz(z, Q, T);

function plotNep(z, p);

[m,n] = size(p);

col={'b','r','k','g','m','c','y'};

% ------------------------------
FIG=mjfig;
AX=linlinax;

f=[];
for i=1:m
	plotstairs(z, p(i,:) +eps, col{i});
end

legendcell={};
for i=1:m
   legendcell={legendcell{:}, ['NEP ',num2str(i)]};
end


legend(legendcell);
hold on
title('siimpl - mean impact parameter')
xlabel('Depth (Å)');
ylabel('NEP (1/Å)');

return
if max(p) == 0, return, end

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

