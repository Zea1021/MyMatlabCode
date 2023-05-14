function [finish_d, totalEner, mt] = calDis(p1_chrom, m1_chrom)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3;

mt = zeros(1, M1);      
finish_d = zeros(1, N); 
totalEner = 0;
for i = 1 : N
    job = p1_chrom(1, i);
    m = m1_chrom(1, i);
    
    startTime = mt(m);
    endTime = startTime + Time1(job, m);
    mt(m) = endTime;  
    finish_d(job) = endTime;           
    totalEner = totalEner + CM1(m) * Time1(job, m); 
end
