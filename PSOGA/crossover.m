function [newChrom1, newChrom2] = crossover(chrom1, chrom2, n, total_L, C, Pc2)
%% operation sequence
os1 = chrom1(1 : total_L);
os2 = chrom2(1 : total_L);
temp1 = os1;
temp2 = os2;
% os cross
L = length(temp1);
position = randperm(L, C);        
for i = 1 : C
    if rand < Pc2
        temp1(position(i)) = os2(position(i));
        temp2(position(i)) = os1(position(i));
    end
end
% os repair
temp1_1 = repair(temp1(1 : n), os1(1 : n));
temp1_2 = repair(temp1(n + 1 : end), os1(n + 1 : end));
newOS1 = [temp1_1, temp1_2];
temp2_1 = repair(temp2(1 : n), os2(1 : n));
temp2_2 = repair(temp2(n + 1 : end), os2(n + 1 : end));
newOS2 = [temp2_1, temp2_2];
%% workstation assign
wa1 = chrom1(total_L + 1 : end);
wa2 = chrom2(total_L + 1 : end);
[newWA1,newWA2]=duodiancross(wa1,wa2);
%% newChrom
newChrom1 = [newOS1, newWA1];
newChrom2 = [newOS2, newWA2];