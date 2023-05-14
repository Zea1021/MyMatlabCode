function newChrom = repair(oldChrom, feasible)
% ��¼
old = oldChrom; 
fea = feasible;

L = length(feasible);
for j = 1 : L
    pos1 = find(fea == old(j), 1);%�ҵ� fea �е�һ������ old��j) ������λ��
    if pos1 > 0
        old(j) = 0;
        fea(pos1) = 0;
    end
end
for j = 1 : L
    if old(j) ~= 0 %����Ļ���
        pos1 = find(oldChrom == old(j), 1);
        pos2 = find(fea, 1);%���� fea �е�һ����Ϊ 0 �Ļ����λ��, �� oldChrom ȱʧ�Ļ���
        oldChrom(pos1) = fea(pos2);%��ȱʧ�Ļ����޲�����Ļ���
        fea(pos2) = 0;
    end
end
newChrom = oldChrom;