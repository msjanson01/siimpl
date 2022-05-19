% function gallralinje(lh);

function out = halfstairs(lh);


if isequal(get(lh,'type'),'line')
   
  		x = get(lh,'xData');
		y = get(lh,'yData');
            
		N = length(x);
      
      N = 4*(floor(N/4));
      x = x(1:N);
      y = y(1:N);
      
		x = sort( [x(1:4:N) x(4:4:N)] );
      
      ym = (y(1:4:N) + y(3:4:N)) / 2;
      
      y = zeros(size(x));
      
      y(1:2:N/2) = ym;
      y(2:2:N/2) = ym;
      
		set(lh,'xData',x);
		set(lh,'yData',y);
     	out=1;
  
  else
      out=-1;
  end
