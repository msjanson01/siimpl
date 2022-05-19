% function Se = Se_siimpl(c, E,  Ef);
% A1 = c(1);
% p  = c(2);
% p2 = c(3);
% A3 = c(4);
% A4 = c(5);
% A5 = c(6);


function Se = Se_siimpl(c, E,  Ef);

A1 = c(1);
p  = c(2);
p2 = c(3);
A3 = c(4);
A4 = c(5);
A5 = c(6);

f = zeros(size(E));
i = find(E > 1.0);
f(i) = exp(-3*(log(Ef)./log(E(i))).^4);

Se_lo  = A1 * E.^( p + f * p2 );

Se_hi = A3./(E/1000) .* log( 1 + A4./(E/1000) + A5*(E/1000) );

Se = (Se_lo.^-1 + Se_hi.^-1).^-1;