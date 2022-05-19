% File name: Se_1H_SiC.m
% Experimental stopping data from M.S. Janson et al. to be published
%
% C and Si stopping are the ones for B and Al, respectively
% Only low-velocity parameters are defined here


% ------- Electronic Stopping Data -------
Efermi 35; % (keV/amu)

%		ion		C		Si
% ---------------------------------------
esp.A1	4.2 	5.5 	4.0;
esp.p   0.372	0.45	0.52;
esp.p2	0   	-0.05	0;

esp.fl	1.0		1.0 	1.0;
esp.s	0.3	0.4	    0.5;

% NOTE. The carbon atoms must earlier have been defined as 
% atom #1, and Silicon as atom #2 (in Z2 instruction)
