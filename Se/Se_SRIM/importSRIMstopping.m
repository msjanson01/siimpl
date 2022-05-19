% function [E, Se, Sn] = importSRIMstopping(Z1, M1, Z2, M2, N, E);

function [E, Se, Sn] = importSRIMstopping(Z1, M1, Z2, M2, N, E);

if (create_SRIM_table(Z1, M1, Z2, M2, N, E));
    [E, Se, Sn] = importSRfile('SR_out.txt');
else
    E = 0;
    Se = 0;
    Sn = 0;
end
