function [indiv1,indiv2]=binarycross(indiv1,indiv2)
L=length(indiv1);
indiv1sub=indiv1;
indiv2sub=indiv2;
binary=randi([0 1],1,L);
index=logical(binary==1);
indiv1(index)=indiv2sub(index);
indiv2(index)=indiv1sub(index);

