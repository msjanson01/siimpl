% fit Se_siimpl to imported SRIM stopping

function [A1, p, p2, A3, A4, A5] = fit_Se_siimpt(E_srim, Se_srim, Ef, A1, p, p2, A3, A4, A5);


% find least square fit of a pearson distribution to the data
options=optimset('diagnostics','off',...
                 'display','off',... % else final
                 'tolfun',1e-4,...
                 'maxiter',500,...
                 'maxfunevals',10000 ...
                ... 'jacobian','off',...
                ... 'LargeScale','off'...
                );
              
          
% start paramaters
c_start = [A1, p, p2, A3, A4, A5];

XDATA = E_srim;
YDATA = Se_srim;

% fit is done on a log versus log scale
[c,resnorm,residual] = lsqcurvefit(inline('log(Se_siimpl(A, B, C))'), c_start, XDATA, log(YDATA),[],[],options, Ef);

A1 = c(1);
p =  c(2);
p2 = c(3);
A3 = c(4);
A4 = c(5);
A5 = c(6);
