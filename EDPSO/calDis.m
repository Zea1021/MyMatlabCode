function [finish_d, totalEner, mt] = calDis(p1_chrom, m1_chrom)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3;

mt = zeros(1, M1);      % completion time of machine
finish_d = zeros(1, N); % completion time of dis
totalEner = 0; % energy
for i = 1 : N
    job = p1_chrom(1, i);
    m = m1_chrom(1, i);
    mt(m) = mt(m) + Time1(job, m);   
    finish_d(job) = mt(m);           
    totalEner = totalEner + CM1(m) * Time1(job, m);
end
