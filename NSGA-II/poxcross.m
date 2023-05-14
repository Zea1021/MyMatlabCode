function [indiv1Sub,indiv2Sub]=poxcross(indiv1,indiv2)
L=length(indiv1);
pos=randperm(L,2);
pos=sort(pos);
S1=indiv1;
S2=indiv2;
Chunk1=indiv1(1,pos(1):pos(2));
Chunk2=indiv2(1,pos(1):pos(2));
for i=1:length(Chunk1)
    pos1=find(S2==Chunk1(i),1);
    S2(pos1)=[];
    pos2=find(S1==Chunk2(i),1);
    S1(pos2)=[];
end
indiv1Sub=[S2(1,1:pos(1)-1),Chunk1,S2(1,pos(1):end)];
indiv2Sub=[S1(1,1:pos(1)-1),Chunk2,S1(1,pos(1):end)];