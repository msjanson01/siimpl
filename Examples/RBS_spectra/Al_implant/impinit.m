% File name: impinit.m

% Crystal, target data
readsimfile D:\siimpl\crystals\4H_SiC.txt

% Electronic stopping
readsimfile D:\siimpl\Se\SiC\Se_Al_SiC.m

% --- Implanted ion data  --------
Z1	    13;
M1   	27; % (amu)
E0	    200; % (keV)

dose    5e14; % (cm-2)

% Implantation direction
impdirrotyrotz	8 0; %  (degrees) 

% --- Damage model ------
damage.cascade 1;

% --- MC statistics ---------
Nions	10000; % number of pseudo ions

layerwidth 4000; % (Å)
Ndistz 60; % division of layer width
