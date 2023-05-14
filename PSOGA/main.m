close all
clear
clc

global n p k MD MP ct Time total_L total_p;
load ../DataSet/4-4.mat
total_L = n + n * sum(k);
total_p= n * sum(k);
%% para
Penum = 120; 
C1 = 2; C2 = 2;
Wmax = 0.9; Wmin = 0.4;
Pc1 = 0.9; C = 40;
Pc2 = 0.6; Pm = 0.1;
Itermax = 500;
%% popInit
Pop=initPop(Penum);
%% calFit
Fit = calFit(Pop);
%% xbest
xbest = Pop;
xbestFit = Fit;
%% gbest
[~, pareto] = nonRank(xbestFit);
gbest = xbest(pareto, :);
gbestFit = xbestFit(pareto,:);
%% speed
V = rand(Penum, 2 * total_L);
%% loop
for i = 1 : Itermax
    %% GA
    subPop = [];
    for j = 1 : 2 : Penum
        if rand < Pc1
            chrom1 = Pop(j, :);
            chrom2 = Pop(j + 1, :);
            idx = randi(size(gbest, 1));
            bestChrom = gbest(idx, :);
            % crossover
            subPop(end+1,:) = crossover(chrom1, bestChrom, n, total_L, C, Pc2);
            subPop(end+1,:) = crossover(chrom2, bestChrom, n, total_L, C, Pc2);
        end
    end
    % mutation
    subPop = aberrance(subPop,Pm);
    subFit = calFit(subPop);
    % merge
    temp_fit = [Fit; subFit];
    temp_pop = [Pop; subPop];
    [temp_fit, idx] = unique(temp_fit, 'rows'); 
    temp_pop = temp_pop(idx, :);
    [Rank, pareto] = nonRank(temp_fit);
    Distance = crowdingDistance(temp_fit, Rank);
    [~, idx] = sortrows([Rank, Distance], [1, -2]);
    Pop = temp_pop(idx(1 : Penum), :);
    Fit = temp_fit(idx(1 : Penum), :);
    gbest = temp_pop(pareto, :);
    gbestFit = temp_fit(pareto, :);
    %% PSO
    W = Wmax - (Wmax - Wmin)/Itermax * i;
    subPop = zeros(Penum, 2 * total_L);
    subFit = zeros(Penum, 2);
    for j = 1 : Penum
        % update V
        V(j, :) = W * V(j, :) + C1 * rand * (xbest(j, :) - Pop(j, :)) + C2 * rand * (gbest(unidrnd(size(gbest, 1)), :) -  Pop(j, :));
        % operation sequence, os
        os = Pop(j, 1 : total_L);
        % new os
        newOS = [sortOfThePosition(os(1 : n), V(j, 1 : n)), sortOfThePosition(os(n + 1 : end), V(j, n + 1 : total_L))];
        % workstation assignment, wa
        wa = Pop(j, total_L + 1 : end);
        % new wa
        newWA = round(wa + V(j, total_L + 1 : end));
        for ii = 1 : n
            if newWA(ii) < 1 || newWA(ii) > MD
                newWA(ii) = unidrnd(MD);
            end
        end
        for ii = n + 1 : total_L
            if newWA(ii) < 1 || newWA(ii) > MP
                newWA(ii) = unidrnd(MP);
            end
        end
        % update new chrom
        subPop(j, :) = [newOS, newWA];
        subFit(j, :) = calFit(subPop(j, :));
        if (isdominate(subFit(j, :), xbestFit(j, :)))
            xbest(j, :) = subPop(j, :);
            xbestFit(j, :) = subFit(j, :);
        end
    end
    % update
    temp_fit = [Fit; subFit];
    temp_pop = [Pop; subPop];
    [temp_fit, idx] = unique(temp_fit, 'rows');
    temp_pop = temp_pop(idx, :);
    [Rank, pareto] = nonRank(temp_fit);
    Distance = crowdingDistance(temp_fit, Rank);
    [~, idx] = sortrows([Rank, Distance], [1, -2]);
    Pop = temp_pop(idx(1 : Penum), :);
    Fit = temp_fit(idx(1 : Penum), :);
    gbest = temp_pop(pareto, :);
    gbestFit = temp_fit(pareto, :);
end