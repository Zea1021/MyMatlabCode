clear
close all
clc

global n p k MD MP ct Time total_L total_p;
load ../DataSet/4-4.mat
%% para
PS=120;
Pc=0.8;
Pm=0.2;
maxgen=500;
total_L = n + n * sum(k); 
total_p = n * sum(k);     
%% init and cal
Pop = initPop(PS);
Fit = calFit(Pop);
%% loop
for gen=1:maxgen
    %% cross
    SubPop = cross(Pop,Pc);
    %% mutation
    SubPop = aberrance(SubPop,Pm);
    %% cal
    SubFit = calFit(SubPop);
    %% update
    BigPop = [Pop; SubPop];
    BigFit = [Fit; SubFit];
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