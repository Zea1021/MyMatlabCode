function [Pop,Fit,E,FitE]=FFOsearch(Pop,Fit,E,FitE,b1,PS)
global n p k MD MP ct Time total_L total_p;
BigPop=[];
BigFit=[];
PopSize = size(Pop,1);
a1 = 0.3; a2 = 7.5;
a = a1 + rand * (a2 - a1);
F = 0.5;
%% Smell-searching
for i=1:PopSize
    chrom = Pop(i,:);
    fit = Fit(i,:);
    
    tempPop=[];
    tempFit=[];
    while size(tempPop,1)<(b1)
        if rand<(F)
            idea=1;
        else
            idea=2;
        end
        switch idea
            case 1
            val = chrom(1,1:n);
            valsub = swapOrInsert(val);
            newChrom = [valsub,chrom(1,n+1:end)];
            case 2
            val = chrom (1,n+1:total_L);
            valsub = swapOrInsert(val);
            newChrom = [chrom(1,1:n),valsub,chrom(1,total_L+1:end)];
        end
        % cal
        newFit = cal(newChrom);
        % judge
        if (newFit(1,1)<fit(1,1))&&(newFit(1,2)<fit(1,2))
            chrom=newChrom;
            fit=newFit;
        elseif (newFit(1,1)>fit(1,1))&&(newFit(1,2)<=fit(1,2))&&rand<exp((-(newFit(1,1)-fit(1,1)))/a)
            chrom=newChrom;
            fit=newFit;
        elseif (newFit(1,1)<=fit(1,1))&&(newFit(1,2)>fit(1,2))&&rand<exp((-(newFit(1,2)-fit(1,2)))/(n*a))
            chrom=newChrom;
            fit=newFit;
        elseif (newFit(1,1)>fit(1,1))&&(newFit(1,2)>fit(1,2))&&((newFit(1,1)-fit(1,1))/(newFit(1,1)))<((newFit(1,2)-fit(1,2))/(newFit(1,2)))&&rand<exp((-(newFit(1,1)-fit(1,1)))/a)
            chrom=newChrom;
            fit=newFit;
        elseif (newFit(1,1)>fit(1,1))&&(newFit(1,2)>fit(1,2))&&((newFit(1,1)-fit(1,1))/(newFit(1,1)))>=((newFit(1,2)-fit(1,2))/(newFit(1,2)))&&rand<exp((-(newFit(1,1)-fit(1,1)))/a)
            chrom=newChrom;
            fit=newFit;
        end
        % record
        tempPop(end+1:end+1,:) = newChrom;
        tempFit(end+1:end+1,:) = newFit;
    end
    % record
    BigPop(end+1:end+size(tempPop,1),:)=tempPop;
    BigFit(end+1:end+size(tempPop,1),:)=tempFit;
end
%% Vision-searching
BigPop(end+1:end+size(E,1),:) = E;
BigFit(end+1:end+size(E,1),:) = FitE;
[BigFit, idx] = unique(BigFit, 'rows'); 
BigPop = BigPop(idx, :);
[Rank, pareto] = nonRank(BigFit);
E = BigPop(pareto, :);
FitE  = BigFit(pareto, :);
to_remove = size(BigPop, 1) - PS; 
if to_remove > 0
    Distance = crowdingDistance(BigFit, Rank);
    [~, idx] = sortrows([Rank, Distance], [1, -2]);
    Pop = BigPop(idx(1 : PS), :);
    Fit = BigFit(idx(1 : PS), :);
end