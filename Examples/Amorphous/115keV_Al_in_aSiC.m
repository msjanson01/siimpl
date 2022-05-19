% File name: 115keV_Al_in_SiC.m

% Crystal, target data
readsimfile D:\siimpl\crystals\4H_SiC.txt

randomtarget

% Electronic stopping
readsimfile D:\siimpl\Se\SiC\Se_Al_SiC.m

% --- Implanted ion data  --------
Z1	    13;
M1   	27; % (amu)
E0	    115; % (keV)

dose    5e14; % (cm-2)

% Implantation direction
impdirrotyrotz	8 0; %  (degrees) 

% --- Damage model ------
damage.cascade 0;

% --- MC statistics ---------
Nions	10000; % number of pseudo ions

layerwidth 5000; % (Å)
Ndistz 50; % division of layer width
