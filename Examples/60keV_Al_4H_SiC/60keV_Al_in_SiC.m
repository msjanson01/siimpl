% File name: 60keV_Al_in_SiC.m

% Crystal, target data
readsimfile D:\siimpl\crystals\4H_SiC.txt

% Electronic stopping
readsimfile D:\siimpl\Se\SiC\Se_Al_SiC.m

% --- Implanted ion data  --------
Z1	    13;
M1   	27; % (amu)
E0	    60; % (keV)

dose    1e14; % (cm-2)

% Implantation direction
impdirrotyrotz	8 0; %  (degrees) 

% --- Damage model ------
damage.cascade 1;

% --- MC statistics ---------
Nions	100000; % number of pseudo ions

layerwidth 10000; % (Å)
Ndistz 100; % division of layer width
