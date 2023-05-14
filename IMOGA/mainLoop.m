function [Elite, FitE] = mainLoop(filename)
clearvars -except filename;
global n p k MD MP ct Time total_L total_p;
prefix = '../DataSet/';
suffix = '.mat';
datafile = [prefix filename suffix]; % filename = '4-4';
load (datafile);
%% para
PS = 120;
Pc = 0.9;
Pm = 0.8;
F = 0.2;
maxgen = 500;
total_L = n + n * sum(k); 
total_p = n * sum(k);     
%% init
Pop = initPopByRule(PS,PS/4);
Fit = calFit(Pop);
%% elite
[~, pareto] = nonRank(Fit);
Elite = Pop(pareto, :);
FitE  = Fit(pareto, :);
%% loop
for gen=1:maxgen
    %% group 1
    [~, index1] = sort(Fit(:, 1));
    %% group 2
    [~, index2] = sort(Fit(:, 2));
    %% cross
    Pop1 = cross(Pop(index1(1 : PS/2),:),Pc);
    Pop2 = cross(Pop(index2(1 : PS/2),:),Pc);
    %% aberrance
    Pop1 = aberrance1(Pop1,Pm,F);
    Pop2 = aberrance2(Pop2,Pm,F);
    %% cal
    Fit1 = calFit(Pop1);
    Fit2 = calFit(Pop2);
    %% update
    SubPop = [Pop1; Pop2];
    SubFit = [Fit1; Fit2];
    [~, pareto] = nonRank(SubFit);
    Elite(end+1:end+length(pareto), :) = SubPop(pareto, :);
    FitE(end+1:end+length(pareto), :)  = SubFit(pareto, :);
    %% local search
    [newPop, newFit] = localSearch(Elite, FitE, F);
    %% next
    BigPop = [Pop; SubPop; newPop];
    BigFit = [Fit; SubFit; newFit];
    [BigFit, idx] = unique(BigFit, 'rows');
    BigPop = BigPop(idx, :);
    [Rank, pareto] = nonRank(BigFit);
    Elite = BigPop(pareto, :);
    FitE  = BigFit(pareto, :);
    to_remove = size(BigPop, 1) - PS; 
    if to_remove > 0
        Distance = crowdingDistance(BigFit, Rank);
        [~, idx] = sortrows([Rank, Distance], [1, -2]);
        Pop = BigPop(idx(1 : PS), :);
        Fit = BigFit(idx(1 : PS), :);
    end
end