function [new_p1_chrom1, new_p2_chrom1, new_p1_chrom2, new_p2_chrom2] = preCross(p1_chrom1, p2_chrom1, p1_chrom2, p2_chrom2)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3; 
temp = randi([0 1], 1, N);
while sum(temp) == 0 || sum(temp) == N
    temp = randi([0 1], 1, N);
end
jobSet = 1 : N;
set1 = jobSet(temp == 0);
set2 = jobSet(temp == 1);
% ≥ı ºªØ
new_p1_chrom1 = zeros(1, N);
new_p2_chrom1 = zeros(1, SH);
new_p1_chrom2 = zeros(1, N);
new_p2_chrom2 = zeros(1, SH);

for i = 1: N
    if (ismember(p1_chrom1(i), set1))
        new_p1_chrom1(i) = p1_chrom1(i);
    end
    if (ismember(p1_chrom2(i), set1))
        new_p1_chrom2(i) = p1_chrom2(i);
    end
end
for i = 1 : N
    if (ismember(p1_chrom2(i), set2))
        idx = find(new_p1_chrom1 == 0, 1);
        new_p1_chrom1(idx) = p1_chrom2(i);
    end
    if (ismember(p1_chrom1(i), set2))
        idx = find(new_p1_chrom2 == 0, 1);
        new_p1_chrom2(idx) = p1_chrom1(i);
    end
end

for i = 1: SH
    if (ismember(floor(p2_chrom1(i)/100), set1))
        new_p2_chrom1(i) = p2_chrom1(i);
    end
    if (ismember(floor(p2_chrom2(i)/100), set1))
        new_p2_chrom2(i) = p2_chrom2(i);
    end
end
for i = 1 : SH
    if (ismember(floor(p2_chrom2(i)/100), set2))
        idx = find(new_p2_chrom1 == 0, 1);
        new_p2_chrom1(idx) = p2_chrom2(i);
    end
    if (ismember(floor(p2_chrom1(i)/100), set2))
        idx = find(new_p2_chrom2 == 0, 1);
        new_p2_chrom2(idx) = p2_chrom1(i);
    end
end