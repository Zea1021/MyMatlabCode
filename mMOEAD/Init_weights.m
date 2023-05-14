function [weights,neighbour] = Init_weights(H,M,T)
N=factorial(H+M-1)/(factorial(M-1)*factorial(H)); 
Popsize=int16(N);  
weights=zeros(Popsize,M);  % init weights
h=(0:H)/H; 
count=1;
for i=1:H+1 
   for j=1:H+1
         if h(i)+h(j)==1 
            temp=[h(i) h(j)]; 
            weights(count,:)=temp;
            count=count+1;
         end
   end
end
distance = zeros(Popsize,Popsize); 
neighbour=zeros(Popsize,T); 
for i=1:Popsize
        for j=i+1:Popsize
            A=weights(i,:);B=weights(j,:);
            distance(i,j)=(A-B)*(A-B)';  
            distance(j,i)=distance(i,j); 
        end
        [~,sindex]=sort(distance(i,:));  
        neighbour(i,:)=sindex(1:T); 
end
end