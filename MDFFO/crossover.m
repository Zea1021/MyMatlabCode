function SubE=crossover(E,FitE,b2,F)
global n p k MD MP ct Time total_L total_p;
Esize = size(E,1);
SubE = [];
if Esize <= 2
    indiv=E(unidrnd(Esize),:);
    for i=1:ceil(n*b2)
        if rand<(F)
            idea=1;
        else
            idea=2;
        end
        switch idea
            case 1
                val=indiv(1,1:n);
                valsub=swapOrInsert(val);
                newChrom=[valsub,indiv(1,n+1:end)];
            case 2
                val=indiv(1,n+1:total_L);
                valsub=swapOrInsert(val);
                newChrom=[indiv(1,1:n),valsub,indiv(1,total_L+1:end)];
        end
        SubE(end+1:end+1,:)=newChrom;
    end
else
    RankE=ones(Esize,1); % init for the same level
    DistE=crowdingDistance(FitE,RankE);
    idx=find(DistE==inf);
    tempDistE=DistE;
    tempDistE(idx)=[];
    DistE(idx)=2*max(tempDistE);
    P=cumsum(DistE)/sum(DistE);
    for i=1:ceil(n*b2)
        index1 = RouletteWheelSelection(P);
        index2 = RouletteWheelSelection(P);
        while index1==index2
            index1 = RouletteWheelSelection(P);
            index2 = RouletteWheelSelection(P);
        end
        val1=E(index1,:);
        val2=E(index2,:);
        
        Subval=toCross(val1,val2);
        SubE(end+1:end+2,:)=Subval;
    end
end