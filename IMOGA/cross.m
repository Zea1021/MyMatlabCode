function SubChrom = cross(Chrom, Pc)
popSize=size(Chrom,1);
index=randperm(popSize);
pop=Chrom(index,:);
SubChrom=[];
for i=1:2:popSize
    if rand<Pc          
        val1=pop(i,:);
        val2=pop(i+1,:);
        Subval=toCross(val1,val2);
        SubChrom(end+1:end+2,:)=Subval;
    end
end
    