% plot net defects

[z, Czion]=import_mj('Czion.txt');

% --- interstitials
A=import_mj('CzI.txt');

if isempty(A),disp('Can not open CzI.txt file!'); return; end
[m,n] = size(A);
zi = A(:,1)';
CzI = A(:,2:n)';

% --- vacancies
A=import_mj('CzV.txt');

if isempty(A),disp('Can not open CzV.txt file!'); return; end
[m,n] = size(A);
zv = A(:,1)';
CzV = A(:,2:n)';
% ----------------------------

mjfig;
linlogax;
title('SiiMPL - net point defect distributions')

plotstairs(z,Czion+eps,'k');

%  --  C ---
IC = CzI(1,:);
VC = CzV(1,:);

IV = (IC - VC);

Inet=IV; 
Vnet=-IV;

Inet(find(Inet<0)) = 0;
ICnet = Inet;

Vnet(find(Vnet<0)) = 0;
VCnet = Vnet;

plotstairs(z,Inet+eps,'r');
plotstairs(z,Vnet+eps,'g');

%  --  Si ---
ISi = CzI(2,:);
VSi = CzV(2,:);

IV=(ISi - VSi);

Inet=IV; 
Vnet=-IV;

Inet(find(Inet<0)) = 0;
ISinet = Inet;

Vnet(find(Vnet<0)) = 0;
VSinet = Vnet;

plotstairs(z,Inet+eps,'b');
plotstairs(z,Vnet+eps,'m');
% --------------

% Excess interstitials
Iex=( (ICnet -VCnet) - (ISinet -VSinet));

ICex = Iex; 
ISiex = -Iex;

ICex(find(ICex<0)) = 0;
ISiex(find(ISiex<0)) = 0;

l=plot(z,ICex+eps,'.r');
set(l,'linestyle','none','color','r');

l=plot(z,ISiex+eps,'.k');
set(l,'linestyle','none','color','k');
% --------------


legend('ion','net C Int.','net C Vac.',...
   'net Si Int.','net Si Vac.',...
    'Excess C','Excess Si');
xlabel('Depth (Å)');
ylabel('Concentration (cm^{-3})');

win=axis;
axis([win(1:2) 1e17 1e21]);