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



























% %测试  function Distance=Genetic_searching(FitnVE)
% XOVR=0.8;%交叉概率
% MUTA=0.1;%变异概率
% popSize=size(E,1);%种群大小
% n=size(E,2);%产品数量
% numObj=size(FitnVE,2);%目标函数个数
% Distance=zeros(popSize,1);%用于存储拥挤距离
% 
% %如果E中只有一个个体
% if size(E,1)<=1
%         E(2,:)=randperm(n);
%         FitnVE(2,:)=calculate_fitness(E(2,:),T,m,s,r);
% end
% 
% %计算E中个体的拥挤距离
% for i=1:numObj
%     data=FitnVE(:,i);
%     [~,index]=sort(data); %对data即适应度值进行升序排序
%     Distance([index(1),index(end)])=inf;%其两端的个体拥挤距离为inf
%     ii=2;
%     while ii<length(index)      
%         Distance(index(ii))=Distance(index(ii))+min(inf,(data(index(ii+1))-data(index(ii-1)))/(max(data)-min(data)));
%         ii=ii+1;
%     end
% end
% 
% %%
% %基于拥挤距离根据轮盘赌规则从E中选出两个父代个体用于进行遗传操作
% %拥挤距离越大，被选中的概率越大
% index=isfinite(Distance);
% if nnz(index)==0    %即Distance中的元素全为inf
%     Distance=ones(size(Distance,1),1);
% else
%     pos=find(Distance==inf);
%     Distancetemp=Distance;
%     Distancetemp(pos)=[];
%     if nnz(Distancetemp)==0     %如果Distance的元素除了是inf就是0
%         Distance(find(Distance==0))=1;
%         Distance(pos)=10;
%     else
%         Distance(pos)=2*max(Distancetemp);%将inf重新赋值为Distance最大值的倍数，这里为两倍
%     end
% end
% sum_Distance=sum(Distance);
% value=Distance./sum_Distance;
% value=cumsum(value);
% temp=[];
% FitnVtemp=[];
% while size(temp,1)<b2*n
%     numpopSub=2;%轮盘赌选择出的个体数
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
%     %此时ESub中含有两个个体，将作为父代进行遗传搜索
%    
%     %%
%     
%     %交叉
%     ESub=across(ESub,XOVR);
%     %变异
%     ESub=aberrance(ESub,MUTA);
%     
%     %%
%     
%     %计算适应度值
%     FitnVESub=calculate_fitness(ESub,T,m,s,r);
%     %生成临时解集
%     temp=[temp;ESub];
%     FitnVtemp=[FitnVtemp;FitnVESub];
% end

