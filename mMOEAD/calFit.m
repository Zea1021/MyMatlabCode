function [fit] = calFit(p1_chrom, m1_chrom, p2_chrom, m2_chrom)
[finish_d, totalEner] = calDis(p1_chrom, m1_chrom);
[finish_p, totalEner] = calPro(p2_chrom, m2_chrom, finish_d, totalEner);
[finish_a, totalEner] = calAss(finish_p, totalEner);
fit = [finish_a totalEner];