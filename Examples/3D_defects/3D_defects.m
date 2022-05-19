% File name: 3D_ion.m

% Crystal, target data
readsimfile D:\siimpl\crystals\4H_SiC.txt

% Electronic stopping
readsimfile D:\siimpl\Se\SiC\Se_Al_SiC.m

% --- Implanted ion data  --------
Z1	    13;
M1   	27; % (amu)
E0	    60; % (keV)

dose    1; % (cm-2);

% Implantation direction
impdirrotyrotz	8 0; %  (degrees)

Rimplant 0 0 0; % (Å)
implantarea 1000 1000; % (Å*Å)


% --- Damage model ------
damage.cascade 1;

% --- MC statistics ---------
Nions	10; % number of pseudo ions

layerwidth 10000; % (Å)
Ndistz 100; % division of layer width

saveIVR3