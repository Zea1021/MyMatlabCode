function newChrom = repair(oldChrom, feasible)
% 记录
old = oldChrom; 
fea = feasible;

L = length(feasible);
for j = 1 : L
    pos1 = find(fea == old(j), 1);%找到 fea 中第一个等于 old（j) 的数的位置
    if pos1 > 0
        old(j) = 0;
        fea(pos1) = 0;
    end
end
for j = 1 : L
    if old(j) ~= 0 %多余的基因
        pos1 = find(oldChrom == old(j), 1);
        pos2 = find(fea, 1);%查找 fea 中第一个不为 0 的基因的位置, 即 oldChrom 缺失的基因
        oldChrom(pos1) = fea(pos2);%用缺失的基因修补多余的基因
        fea(pos2) = 0;
    end
end
newChrom = oldChrom;