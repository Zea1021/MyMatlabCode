function Distance=crowdingDistance(FitnV,Rank)
popSize=size(FitnV,1);
numObj=size(FitnV,2);
Distance=zeros(popSize,1);
num=size(unique(Rank),1);
for i=1:num % for each front
    front=find(Rank==i);
    Distancetemp=zeros(size(front,1),1);
    FitnVfront=FitnV(front,:);
    for j=1:numObj                   % for each obj
        data=FitnVfront(:,j);        
        [~,index]=sort(data);        % sort
        Distancetemp([index(1),index(end)])=inf; % set edge value with inf
        ii=2;
        while ii<length(index)     
            Distancetemp(index(ii))=Distancetemp(index(ii))+min(inf,(data(index(ii+1))-data(index(ii-1)))/(max(data)-min(data)));
            ii=ii+1;
        end
    end
    Distance(front)=Distancetemp;
end