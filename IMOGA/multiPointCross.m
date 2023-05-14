function [indivsub1,indivsub2]=multiPointCross(indiv1,indiv2)
L=length(indiv1);
indivsub1=indiv1;
indivsub2=indiv2;
pos=randi([0 1],1,L);
pos=find(pos==1);
indivsub1(pos)=indiv2(pos);
indivsub2(pos)=indiv1(pos);