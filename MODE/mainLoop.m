function [E, FitE] = mainLoop(filename)
clearvars -except filename;
global n p k MD MP ct Time total_L total_p;
prefix = '../DataSet/';
suffix = '.mat';
datafile = [prefix filename suffix]; % filename = '4-4';
load (datafile);
%% para
PS=120;
Cr=0.8;
F=0.5;
maxIter=500;
total_L=n+n*sum(k); 
total_p=n*sum(k);   
%% init
Pop = initPop(PS);
Fit = calFit(Pop);
%% elite
[~, pareto] = nonRank(Fit);
E = Pop(pareto, :);
FitE  = Fit(pareto, :);
%% loop
for iter = 1:maxIter
    % aberrance
    SubPop = aberrance(Pop);
    % cross
    SubPop = DEcross(SubPop,E,Cr,F);
    % cal
    SubFit = calFit(SubPop);
    % next
    BigPop = [Pop; SubPop];
    BigFit = [Fit; SubFit];
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
end