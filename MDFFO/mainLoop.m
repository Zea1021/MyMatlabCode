function [E, FitE] = mainLoop(filename)
clearvars -except filename;
global n p k MD MP ct Time total_L total_p;
prefix = '../DataSet/';
suffix = '.mat';
datafile = [prefix filename suffix]; % filename = '4-4';
load (datafile);
%% para
PS=120;
b1=1; b2=2;
maxIter=500;
total_L=n+n*sum(k);
total_p=n*sum(k);  
%% init and cal fit
Pop = initPop(PS);
Fit = calFit(Pop);
%% elite
[~, pareto] = nonRank(Fit);
E = Pop(pareto,:);
FitE = Fit(pareto,:);
%% loop
for iter = 1 : maxIter
    [Pop,Fit,E,FitE] = FFOsearch(Pop,Fit,E,FitE,b1,PS);
    [E,FitE]=GeneticSearch(E,FitE,b2);
end