% moments

function mom = moments(XDATA,YDATA);

   
% calculate moments
xu1=sum(XDATA.*YDATA)/sum(YDATA);
xu2=sum(XDATA.^2.*YDATA)/sum(YDATA);

mu = xu1;
sigma = sqrt(xu2-xu1^2);

if sigma==0
   mom=[mu sigma 0 0];
   return
end

gama = (sum((XDATA-mu).^3.*YDATA)/sum(YDATA))/sigma^3;
beta = (sum((XDATA-mu).^4.*YDATA)/sum(YDATA))/sigma^4;

mom = [mu sigma gama beta];

