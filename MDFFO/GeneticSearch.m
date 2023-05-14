function [E, FitE] = GeneticSearch(E, FitE, b2)
SubE = crossover(E,FitE,b2,0.5);
SubE = aberrance(SubE,0.5);
SubFitE = calFit(SubE);
%% update
BigPop = [E; SubE];
BigFit = [FitE; SubFitE];
[BigFit, idx] = unique(BigFit, 'rows'); 
BigPop = BigPop(idx, :);
[~, pareto] = nonRank(BigFit);
E = BigPop(pareto, :);
FitE  = BigFit(pareto, :);



























% %����  function Distance=Genetic_searching(FitnVE)
% XOVR=0.8;%�������
% MUTA=0.1;%�������
% popSize=size(E,1);%��Ⱥ��С
% n=size(E,2);%��Ʒ����
% numObj=size(FitnVE,2);%Ŀ�꺯������
% Distance=zeros(popSize,1);%���ڴ洢ӵ������
% 
% %���E��ֻ��һ������
% if size(E,1)<=1
%         E(2,:)=randperm(n);
%         FitnVE(2,:)=calculate_fitness(E(2,:),T,m,s,r);
% end
% 
% %����E�и����ӵ������
% for i=1:numObj
%     data=FitnVE(:,i);
%     [~,index]=sort(data); %��data����Ӧ��ֵ������������
%     Distance([index(1),index(end)])=inf;%�����˵ĸ���ӵ������Ϊinf
%     ii=2;
%     while ii<length(index)      
%         Distance(index(ii))=Distance(index(ii))+min(inf,(data(index(ii+1))-data(index(ii-1)))/(max(data)-min(data)));
%         ii=ii+1;
%     end
% end
% 
% %%
% %����ӵ������������̶Ĺ����E��ѡ�����������������ڽ����Ŵ�����
% %ӵ������Խ�󣬱�ѡ�еĸ���Խ��
% index=isfinite(Distance);
% if nnz(index)==0    %��Distance�е�Ԫ��ȫΪinf
%     Distance=ones(size(Distance,1),1);
% else
%     pos=find(Distance==inf);
%     Distancetemp=Distance;
%     Distancetemp(pos)=[];
%     if nnz(Distancetemp)==0     %���Distance��Ԫ�س�����inf����0
%         Distance(find(Distance==0))=1;
%         Distance(pos)=10;
%     else
%         Distance(pos)=2*max(Distancetemp);%��inf���¸�ֵΪDistance���ֵ�ı���������Ϊ����
%     end
% end
% sum_Distance=sum(Distance);
% value=Distance./sum_Distance;
% value=cumsum(value);
% temp=[];
% FitnVtemp=[];
% while size(temp,1)<b2*n
%     numpopSub=2;%���̶�ѡ����ĸ�����
% %         pos=randperm(size(E,1),2);
% %         ESub(1,:)=E(pos(1),:);
% %         ESub(2,:)=E(pos(2),:);
%     if size(E,1)==2
%         ESub=E;
%     else
%         for i=1:numpopSub
%             select=find(value>=rand);
%             to_select=E(select(1),:);
%             ESub(i,:)=to_select;
%             FitnVESub(i,:)=FitnVE(select(1),:);
%         end
%     end
%     %��ʱESub�к����������壬����Ϊ���������Ŵ�����
%    
%     %%
%     
%     %����
%     ESub=across(ESub,XOVR);
%     %����
%     ESub=aberrance(ESub,MUTA);
%     
%     %%
%     
%     %������Ӧ��ֵ
%     FitnVESub=calculate_fitness(ESub,T,m,s,r);
%     %������ʱ�⼯
%     temp=[temp;ESub];
%     FitnVtemp=[FitnVtemp;FitnVESub];
% end

