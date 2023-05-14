function [c_p1_chrom, c_m1_chrom, c_p2_chrom, c_m2_chrom] = evolution(p1_chrom1, m1_chrom1, p2_chrom1, m2_chrom1, p1_chrom2, m1_chrom2, p2_chrom2, m2_chrom2)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3;
global PS Pc Pm;
%% cross
[c_p1_chrom] = tmxCross(p1_chrom1, p1_chrom2);
[c_p2_chrom] = tmxCross(p2_chrom1, p2_chrom2);
[c_m1_chrom] = sbxCross(m1_chrom1, m1_chrom2);
[c_m2_chrom] = sbxCross(m2_chrom1, m2_chrom2);
%% mutation
if rand < Pm 
    % dis
    pos1 = randi(N, 1, 1);
    pos2 = randi(N, 1, 1); 
    while (pos1 == pos2)
        pos2 = randi(N, 1, 1); 
    end
    c_p1_chrom([pos1, pos2]) = c_p1_chrom([pos2, pos1]);
    c_m1_chrom(pos1) = randi(M1);
    % repro
    pos1 = randi(SH, 1, 1); 
    pos2 = randi(SH, 1, 1); 
    while (pos1 == pos2)
        pos2 = randi(SH, 1, 1); 
    end
    c_p2_chrom([pos1, pos2]) = c_p2_chrom([pos2, pos1]);
    c_m2_chrom(pos1) = randi(M2);
end