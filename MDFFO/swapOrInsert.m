function indivsub=swapOrInsert(indiv)
L=size(indiv,2);
if rand<0.5
    indivsub=indiv;
    pos=randperm(L,2);
    temp=indiv(pos(1));
    indivsub(pos(1))=indiv(pos(end));
    indivsub(pos(end))=temp;
else
    pos=randperm(L,2);
    Chunk1=indiv(1:min(pos)-1);
    gene=indiv(max(pos));
    indiv(max(pos))=[];
    Chunk2=indiv(min(pos):end);
    indivsub=[Chunk1 gene Chunk2];
end  