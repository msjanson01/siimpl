% File name: 60keV_Al_0001.m

% Crystal, target data
readsimfile D:\siimpl\crystals\4H_SiC.txt

randsurflayer 8;  % (Å)

% Electronic stopping
readsimfile D:\siimpl\Se\SiC\Se_Al_SiC.m

% --- Implanted ion data  --------
Z1	    13;
M1   	27; % (amu)
E0	    60; % (keV)

dose    1e14; % (cm-2)

% Implantation direction
impdirUC 0 0 1;  

% --- Damage model ------
damage.cascade 1;
damage.c_a  1.2;

% --- MC statistics ---------
Nions	80000; % number of pseudo ions

layerwidth 10000; % (Å)
Ndistz 100; % division of layer width

saveatdose  4e12  1e13 3e13;