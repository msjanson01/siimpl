% plotIVzFile.m

%  interstitials
A=import_mj('CzI.txt');

if isempty(A),disp('Can not open CzI.txt file!'); return; end
[m,n] = size(A);
zi = A(:,1)';
CzI = A(:,2:n)';

% vacancies
A=import_mj('CzV.txt');

if isempty(A),disp('Can not open CzV.txt file!'); return; end
[m,n] = size(A);
zv = A(:,1)';
CzV = A(:,2:n)';
% --------------


if ( prod( double(zi==zv)) ~= 1)
   disp('Iz.txt and Vz.txt files are not compatible!');
   return;
end

% ------------------------------
plotCzIV(zi, CzI, CzV);