function [p1_chrom, m1_chrom, p2_chrom, m2_chrom] = initPop(PS)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3;

p1_chrom=zeros(PS,N);% operation
p2_chrom=zeros(PS,sum(sum(H)));
m1_chrom=zeros(PS,N);% workstation
m2_chrom=zeros(PS,sum(sum(H)));

original = [];
for j = 1 : N
    for k = 1 : C
        for l = 1 : H(j, k)
            original = [original j * 100 + k];
        end
    end
end

for i = 1 : PS
    p1_chrom(i, :) = randperm(N, N);
    m1_chrom(i, :) = randi(M1, 1, N);
    
    p2_chrom(i, :) = original(randperm(SH)); 
    m2_chrom(i, :) = randi(M2, 1, SH);
end