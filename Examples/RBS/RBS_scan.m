% File name: RBS_scan.m
rotz = 0;
scan = [0:0.1:1. 1.2:0.2:4];

RBSnormal = 0.270;

% NEP integration interval
a = 250;  % (A)
b = 950;  % (A)

yield = [];

for i = 1:length(scan);
  	% ---- create siimpl init file;
   
   	fid = fopen('init.txt','w');
   
    fprintf(fid,'readsimfile d:\\siimpl\\examples\\RBS\\2MeV_He_rbs.m\r\n');
      
    roty = scan(i);
    fprintf(fid,'impdirrotyrotz %e %e;\r\n', [roty rotz]);

    fclose(fid);
   
   	% --- start siimpl simulation
   	tic
	!d:\siimpl\siimpl init.txt  
   	time=toc;
   	disp(['Simulation time: ',num2str(toc,2)]);disp('');

    [z, NEP1, NEP2] = import_mj('NEP.txt');
    
    % Integrate I Si NEP (NEP2) from a to b
    dz = z(2) - z(1);
    j = find( (z >= a) & (z <= b));
    I = sum(NEP2(j)*dz);
   
    yield = [yield, I /( RBSnormal*(b - a))];
 end
 
 
 save('RBSscan.mat','scan','yield');
 
 axes;
 xlabel('scan (degrees)');
 ylabel('Normalized yield');
 title('4H-SiC, Si-side: from <0001> in (10-10) plane');
 
 plot(scan, yield, 'k-');