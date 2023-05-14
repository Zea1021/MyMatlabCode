function [newPop,newFit]=localSearch(Pop, Fit, F)
Size = size(Pop,1);
newPop = [];
newFit = [];
for i = 1 : Size
    chrom = Pop(i, :);
    fit = Fit(i, :);
    %% PSC
    newChrom1 = PSC(chrom);
    fit1      = cal(newChrom1);
    if isdominate(fit1, fit)
        chrom = newChrom1;
    end
    %% BWL
    newChrom2 = BWL(chrom, F);
    fit2      = cal(newChrom2);
    %% new
    newPop(end+1:end+2, :) = [newChrom1; newChrom2];
    newFit(end+1:end+2, :) = [fit1; fit2];
end