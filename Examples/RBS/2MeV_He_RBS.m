% File name: 2MeV_He_RBS.m

% Crystal, target data
readsimfile D:\siimpl\crystals\4H_SiC.txt

layerwidth 1000; % (Å)

% --- Implanted ion data  --------
Z1	    2;
M1 	    4; % (amu)
E0	    2e3; % (keV)

dose    1e0; % (cm-2)

% First rot around y-axis, the rot around z-axis
impdirrotyrotz 8 8; % (degrees) 

% ------- Electronic Stopping Data -------
% 
%     		ion 	atom1	 atom2
% ----------------------------
esp.A1		1.7;
esp.A2		0.602;
esp.A3		50.5;
esp.A4		2.68;
esp.A5		2.1;

esp.fl		1.0		1.0	    1.0;
esp.s		0.4		0.3	    0.5;


% --- MC statistics ---------
Nions	1000; % number of pseudo ions

Ndistz 40; % division of layer width
