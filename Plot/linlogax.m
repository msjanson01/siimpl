
function varargout=linlogax(varargin);

a=axes;
hold on;

xScale='Linear'; yScale='log';

Title='';
xLabel='';
yLabel='';


set(a,'Box','on','FontSize',16,'LineWidth',1,'MinorGridLineStyle','none');
set(a,'units','norm','position',[  0.100    0.10    0.84    0.82]);
grid off

th=title(Title);xh=xlabel(xLabel);yh=ylabel(yLabel);
set(th,'fontsize',16,'FontWeight','bold');

set(xh,'fontsize',16,'FontWeight','norm');
set(yh,'fontsize',16,'FontWeight','norm');

set(a,'XScale',xScale,'YScale',yScale);

if nargout==1, varargout{1}=a;end