function Chrom=aberrance(Chrom,Pm)
global n p k MD MP ct Time total_L total_p;
num_pop=size(Chrom,1);
for i=1:num_pop
    if rand<Pm           
        val=Chrom(i,:);
        S1=val(1,1:total_L);
        M1=val(1,total_L+1:end);
        
        position=randperm(n,2);
        Chunk_SD=S1(1,1:n);
        Chunk_SD=insert_b(Chunk_SD,position);
        S1=[Chunk_SD,S1(1,n+1:end)];
        
        position=randperm(total_p,2);
        Chunk_SP=S1(1,n+1:total_L);
        Chunk_SP=insert_b(Chunk_SP,position);
        S1=[S1(1,1:n),Chunk_SP];
        Subval=[S1,M1];
        
        Chrom(i,:) = Subval;
    end
end