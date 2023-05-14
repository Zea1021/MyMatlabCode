function SubChrom = cross(Chrom1,Chrom2,F)
global n p k MD MP ct Time total_L total_p;
%% operation
P1_1 = Chrom1(1:n);
P1_2 = Chrom2(1:n);
P2_1 = Chrom1(n+1:total_L);
P2_2 = Chrom2(n+1:total_L);
num1 = ceil(F * n);
num2 = ceil(F * total_p);
pos1 = randperm(n, num1);
pos2 = randperm(total_p, num2);
P1_1(pos1) = P1_2(pos1);
P1_1 = repair(P1_1, Chrom1(1:n));
P2_1(pos2) = P2_2(pos2);
P2_1 = repair(P2_1, Chrom1(n+1:total_L));
%% machine
M1 = Chrom1(total_L+1:end);
M2 = Chrom2(total_L+1:end);
pos = randi([0 1],1,total_L);
M1(pos == 1) = M2(pos == 1);

SubChrom = [P1_1 P2_1 M1];