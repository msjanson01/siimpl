% function create_SRIM_table(Z1, M1,, Z2, M2, N, E);

function out = create_SRIM_table(Z1, M1, Z2, M2, N, E);

if (length(Z2) ~= length(M2)) | (length(Z2) ~= length(N))
    disp('Z2, M2, and N has to be of the same length!');
    out = 0;
    return;
end


% Create SR.IN data file for SRIM S_Module

fid = fopen('SR.IN','w');
   
fprintf(fid,'---Input data for SRIM SR_module created by "create_SRIM_table.m"\r\n');
fprintf(fid,'---Output File Name\r\n');
fprintf(fid,'"SR_out.txt"\r\n');

fprintf(fid,'---Ion data: Z1  M1\r\n');
fprintf(fid,'%d  %f\r\n',[Z1 M1]);

fprintf(fid,'---Target Data: (Solid=0,Gas=1), Density(g/cm3) (not used), Compound Corr.\r\n');
if length(Z2 == 1)
    fprintf(fid,'0    1.0    1.0\r\n');
else
    % NOTE currently no compund correction used, i.e. CC = 1.0!
    fprintf(fid,'0    1.0    1.0\r\n');
end

fprintf(fid,'---Number of target elements\r\n');
fprintf(fid,'%d\r\n',length(Z2));

fprintf(fid,'---Target Elements: (Z), Target name, Stoich, Target Mass(u)\r\n');
for i=1:length(Z2)
    fprintf(fid,'%d "name"  %d   %f\r\n',[Z2(i) N(i), M2(i)]);
end

fprintf(fid,'---Output Stopping Units (1-8): 7 = [eV / (1e15 atoms/cm2)]\r\n');
fprintf(fid,'7\r\n');


fprintf(fid,'---Ion Energy : E-Min(keV), E-Max(keV))\r\n');
fprintf(fid,'0  0\r\n');

for i=1:length(E)
    fprintf(fid,'%d\r\n',E(i));
end

fprintf(fid,'0\r\n');

fclose(fid);


% Execute SR_module
! SRModule.exe

% substitute ',' to '.' in SR_out.txt
! subst_comma SR_out.txt

out = 1;