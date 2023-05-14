function [p1_chrom, p2_chrom, m1_chrom, m2_chrom, fit] = initPopByRule(PS,ruleNum)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3; % ���ȫ�ֱ������⴫��Ƶ������ʱ������
p1_chrom=zeros(PS,N); %��ж������
p2_chrom=zeros(PS,SH);%�ӹ�������
m1_chrom=zeros(PS,N); %��ж������
m2_chrom=zeros(PS,SH);%�ӹ�������
fit = zeros(PS, 2);   % �洢������ʱ����ܵ��ܺ�
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