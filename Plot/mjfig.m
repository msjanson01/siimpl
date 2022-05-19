
function varargout=mjFig;

f=figure;
set(f,'units','norm','position',[0.05    0.05    0.7882    0.8]);
set(f,'PaperUnits','centimeter','PaperOrientation','landscape');
set(f,'PaperType','a4letter');
set(f,'PaperPosition',[3 2 24 16.5]);

if nargout==1, varargout{1}=f;end