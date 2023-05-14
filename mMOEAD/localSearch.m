function [p1_chrom, m1_chrom, p2_chrom, m2_chrom] = localSearch(p1_chrom, m1_chrom, p2_chrom, m2_chrom)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3; 
%% disassembly
pos = randi(N);
job = p1_chrom(pos);
m = 1;
for j = 2 : M1
    if (Time1(job, j) < Time1(job, m))
        m = j;
    end
end
m1_chrom(pos) = m;
%% reprocessing
pos = randi(SH);
signal = p2_chrom(pos);
job = floor(signal/100); 
compo = signal - job * 100; 
idx = find(p2_chrom == signal);
k = find(idx == pos);
m = 1;
for j = 2 : M2
    if (Time2(job, compo, k, j) < Time2(job, compo, k, m))
        m = j;
    end
end
m2_chrom(pos) = m;