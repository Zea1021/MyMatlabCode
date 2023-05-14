function [Chrom,Allot]=initPop(Nind)
global n p k MD MP ct Time total_L total_p;
Chrom=zeros(Nind,2*total_L);            
 
Allot=[];
for ii=1:n
    for iii=1:p
        comp=ii*100+iii;
        Allot(1,end+1:end+k(iii))=comp*ones(1,k(iii));
    end
end

for i=1:Nind
    S_d=randperm(n);
    S_p=Allot(randperm(total_p));
    M_d=randi(MD,1,n);
    M_p=randi(MP,1,total_p);
    
    S=[S_d,S_p,M_d,M_p];
    Chrom(i,:)=S;
end
    
    
        