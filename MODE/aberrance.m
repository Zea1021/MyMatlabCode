function SubPop=aberrance(Pop)
global n p k MD MP ct Time total_L total_p;
PopSize = size(Pop,1);
SubPop  = zeros(PopSize, 2*total_L);
for i = 1 : PopSize
    chrom = Pop(i,:);
    S1 = chrom(1,1:total_L);
    M1 = chrom(1,total_L+1:end);
    
    position=randperm(n,2);
    Chunk_SD=S1(1,1:n);
    Chunk_SD=insert_b(Chunk_SD,position);
    S1=[Chunk_SD,S1(1,n+1:end)];
    
    position=randperm(total_p,2);
    Chunk_SP=S1(1,n+1:total_L);
    Chunk_SP=insert_b(Chunk_SP,position);
    S1=[S1(1,1:n),Chunk_SP];
    SubChrom=[S1,M1];
    
    SubPop(i,:)=SubChrom;
end