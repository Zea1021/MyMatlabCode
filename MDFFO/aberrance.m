function Pop = aberrance(Pop,F)
global n p k MD MP ct Time total_L total_p;
PopSize=size(Pop,1);
for i=1:PopSize
    indiv=Pop(i,:);
    if rand<(F)
        idea=1;
    else
        idea=2;
    end
    switch idea
        case 1
            val=indiv(1,1:n);
            valsub=swapOrInsert(val);
            newchrom=[valsub,indiv(1,n+1:end)];
        case 2
            val=indiv(1,n+1:total_L);
            valsub=swapOrInsert(val);
            newchrom=[indiv(1,1:n),valsub,indiv(1,total_L+1:end)];
    end
    Pop(i,:)=newchrom;
end

            
        
        
