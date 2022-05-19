% function plotStairs(xData, yData);

function varargout=plotStairs(xData, yData, varargin);

if (nargin == 4), 
   AX = varargin{2};
   linecol= varargin{1};
elseif (nargin == 3), 
   linecol= varargin{1};
   AX = gca;
elseif (nargin == 2), 
   AX = gca;
   linecol= 'b';
else 
   disp('Wrong nuber of input statements!');
   return
end

CURRENTFIG=gcf;

FIGHNDL=get(AX,'parent');
set(0,'currentFig',FIGHNDL);
set(FIGHNDL,'CurrentAxes',AX);

[m,n]=size(xData);
if m>n
   xData=xData';
   yData=yData';
end


dx=(xData(2)-xData(1));
N = length(xData);
 
slask = gca;


lineh=stairs([xData, (xData(N) + dx)] - dx/2,...
				 [yData,  yData(N)], linecol);

set(0,'currentFig',CURRENTFIG);

if nargout==1, varargout{1}=lineh; end