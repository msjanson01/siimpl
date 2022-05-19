% plotIVzFile.m

A=import_mj('NEP.txt');

if isempty(A),disp('Can not open Fnep.txt file!'); return; end
[m,n] = size(A);
zi = A(:,1)';
p = A(:,2:n)';

% --------------


% ------------------------------
plotNep(zi, p);