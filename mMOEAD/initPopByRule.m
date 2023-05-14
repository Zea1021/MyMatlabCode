function [p1_chrom, p2_chrom, m1_chrom, m2_chrom, fit] = initPopByRule(PS,ruleNum)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3; % 变成全局变量避免传参频繁减少时间消耗
p1_chrom=zeros(PS,N); %拆卸工序码
p2_chrom=zeros(PS,SH);%加工工序码
m1_chrom=zeros(PS,N); %拆卸机器码
m2_chrom=zeros(PS,SH);%加工机器码
fit = zeros(PS, 2);   % 存储最大完成时间和总的能耗
for i = 1 : ruleNum
    [p1_chrom(i, :), m1_chrom(i, :), p2_chrom(i, :), m2_chrom(i, :), fit(i, :)] = initMethod1(0.5);
end
for i = ruleNum + 1 : 2 * ruleNum
    [p1_chrom(i, :), m1_chrom(i, :), p2_chrom(i, :), m2_chrom(i, :), fit(i, :)] = initMethod2(0.5);
end
for i = 2 * ruleNum + 1 : 3 * ruleNum
    [p1_chrom(i, :), m1_chrom(i, :), p2_chrom(i, :), m2_chrom(i, :), fit(i, :)] = initMethod3(0.5);
end
for i = 3 * ruleNum + 1 : PS
    [p1_chrom(i, :), m1_chrom(i, :), p2_chrom(i, :), m2_chrom(i, :)] = initPop(1);
    fit(i, :) = calFit(p1_chrom(i, :), m1_chrom(i, :), p2_chrom(i, :), m2_chrom(i, :));
end