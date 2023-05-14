function SubPop = DEcross(Pop,Elite,Cr,F)
PopSize = size(Pop,1);
SubPop = [];
for i = 1:PopSize
    if rand<Cr
        chrom1 = Pop(i,:);
        j1 = randi(PopSize);
        while j1 == i
            j1 = randi(PopSize);
        end
        chrom2 = Pop(j1,:);
        j2 = randi(size(Elite, 1));
        chrom3 = Elite(j2,:);
        
        SubChrom1 = cross(chrom1,chrom2,F);
        SubChrom2 = cross(SubChrom1,chrom3,F);
        SubPop(end+1:end+2, :) = [SubChrom1; SubChrom2];
    end
end