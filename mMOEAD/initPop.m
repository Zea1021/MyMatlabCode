function [p1_chrom, m1_chrom, p2_chrom, m2_chrom, fit] = initPop(PS)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3;%变成全局变量避免传参频繁减少时间消耗

p1_chrom=zeros(PS,N);%工序码
p2_chrom=zeros(PS,sum(sum(H)));%工序码
m1_chrom=zeros(PS,N);%机器码
m2_chrom=zeros(PS,sum(sum(H)));%机器码

original = zeros(1, SH);
count = 1;
for j = 1 : N
    for k = 1 : C
        for l = 1 : H(j, k)
            original(count) = j * 100 + k;
            count = count + 1;
        end
    end
end

for i = 1 : PS
    p1_chrom(i, :) = randperm(N, N);
    m1_chrom(i, :) = randi(M1, 1, N);
    
    p2_chrom(i, :) = original(randperm(SH)); 
    m2_chrom(i, :) = randi(M2, 1, SH);
end
%% cal fit
fit = zeros(PS, 2);% 存储最大完成时间和总的能耗
for i = 1 : PS
    fit(i, :) = calFit(p1_chrom(i, :), m1_chrom(i, :), p2_chrom(i, :), m2_chrom(i, :));
end