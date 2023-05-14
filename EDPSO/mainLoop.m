function [elite, elite_fit] = mainLoop(filename)
clearvars -except filename;
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3; 
global PS Pp Pg Pt;
global pbest_p1_chrom pbest_m1_chrom pbest_p2_chrom pbest_m2_chrom pbestFit;
global gbest_p1_chrom gbest_m1_chrom gbest_p2_chrom gbest_m2_chrom gbestFit;
global tbest_p1_chrom tbest_m1_chrom tbest_p2_chrom tbest_m2_chrom tbestFit;
global cbest_p1_chrom cbest_m1_chrom cbest_p2_chrom cbest_m2_chrom cbestFit;
%% loadData
prefix = 'DataSet/';
suffix = '.mat';
datafile = [prefix filename suffix]; % filename = '4-4';
load (datafile);
%% parameter
PS = 120;
Pp = 0.3;
Pg = 0.3;
Pt = 0.2; 
Itermax = 500;
%% initPop
[p1_chrom, m1_chrom, p2_chrom, m2_chrom] = initPop(PS);
%% calFit
fit = zeros(PS, 2);
for i = 1 : PS
    fit(i, :) = calFit(p1_chrom(i, :), m1_chrom(i, :), p2_chrom(i, :), m2_chrom(i, :));
end
%% pbest
pbest_p1_chrom = p1_chrom;
pbest_m1_chrom = m1_chrom;
pbest_p2_chrom = p2_chrom;
pbest_m2_chrom = m2_chrom;
pbestFit       = fit;
%% gbest
[~, pareto]    = nonRank(pbestFit);
gbest_p1_chrom = pbest_p1_chrom(pareto, :);
gbest_m1_chrom = pbest_m1_chrom(pareto, :);
gbest_p2_chrom = pbest_p2_chrom(pareto, :);
gbest_m2_chrom = pbest_m2_chrom(pareto, :);
gbestFit       = pbestFit(pareto, :);
preNum         = size(gbestFit, 2);
%% tbest
[~, idx]       = max(pbestFit(:, 1));
tbest_p1_chrom = pbest_p1_chrom(idx, :);
tbest_m1_chrom = pbest_m1_chrom(idx, :);
tbest_p2_chrom = pbest_p2_chrom(idx, :);
tbest_m2_chrom = pbest_m2_chrom(idx, :);
tbestFit       = pbestFit(idx, :);
%% cbest
[~, idx]       = max(pbestFit(:, 2));
cbest_p1_chrom = pbest_p1_chrom(idx, :);
cbest_m1_chrom = pbest_m1_chrom(idx, :);
cbest_p2_chrom = pbest_p2_chrom(idx, :);
cbest_m2_chrom = pbest_m2_chrom(idx, :);
cbestFit       = pbestFit(idx, :);
%% main loop
for iter = 1 : Itermax
    %% posUpdate
    [c_p1_chrom, c_m1_chrom, c_p2_chrom, c_m2_chrom, c_fit] = posUpdate(p1_chrom, m1_chrom, p2_chrom, m2_chrom);
    %% gbest
    t_p1_chrom = [gbest_p1_chrom; c_p1_chrom];
    t_m1_chrom = [gbest_m1_chrom; c_m1_chrom];
    t_p2_chrom = [gbest_p2_chrom; c_p2_chrom];
    t_m2_chrom = [gbest_m2_chrom; c_m2_chrom];
    t_fit      = [gbestFit; c_fit];
    [t_fit, idx] = unique(t_fit, 'rows'); 
    t_p1_chrom = t_p1_chrom(idx, :);
    t_m1_chrom = t_m1_chrom(idx, :);
    t_p2_chrom = t_p2_chrom(idx, :);
    t_m2_chrom = t_m2_chrom(idx, :);
    [Rank, pareto] = nonRank(t_fit);
    gbest_p1_chrom = t_p1_chrom(pareto, :);
    gbest_m1_chrom = t_m1_chrom(pareto, :);
    gbest_p2_chrom = t_p2_chrom(pareto, :);
    gbest_m2_chrom = t_m2_chrom(pareto, :);
    gbestFit = t_fit(pareto, :);
    %% new population
    [~, idx] = sort(Rank);
    num10 = ceil(PS/10);
    maintain = zeros(1, PS);
    maintain(1 : num10) = idx(1 : num10);
    idx(1 : num10) = [];
    temp = randi(length(idx), 1, PS - num10);
    maintain(num10 + 1 : PS) = idx(temp);
    p1_chrom = t_p1_chrom(maintain, :);
    m1_chrom = t_m1_chrom(maintain, :);
    p2_chrom = t_p2_chrom(maintain, :);
    m2_chrom = t_m2_chrom(maintain, :);
    fit = t_fit(maintain, :);
    %% local search
    rand10 = randperm(PS, ceil(PS/10));
    for i = 1 : length(rand10)
        [n_p1_chrom, n_m1_chrom, n_p2_chrom, n_m2_chrom, n_fit] = localSearch(p1_chrom(rand10(i), :), m1_chrom(rand10(i), :), p2_chrom(rand10(i), :), m2_chrom(rand10(i), :), fit(rand10(i), :));
        % update chrom
        [~, pareto] = nonRank(n_fit);
        idx = randi(length(pareto));
        p1_chrom(rand10(i), :) = n_p1_chrom(pareto(idx), :);
        m1_chrom(rand10(i), :) = n_m1_chrom(pareto(idx), :);
        p2_chrom(rand10(i), :) = n_p2_chrom(pareto(idx), :);
        m2_chrom(rand10(i), :) = n_m2_chrom(pareto(idx), :);
        fit(rand10(i), :)      = n_fit(pareto(idx), :);
        % pbest
        for j = randperm(length(pareto))
            flag = isdominate(n_fit(pareto(j), :), pbestFit(rand10(i), :));
            if flag == 1
                pbest_p1_chrom(rand10(i), :) = n_p1_chrom(pareto(j), :);
                pbest_m1_chrom(rand10(i), :) = n_m1_chrom(pareto(j), :);
                pbest_p2_chrom(rand10(i), :) = n_p2_chrom(pareto(j), :);
                pbest_m2_chrom(rand10(i), :) = n_m2_chrom(pareto(j), :);
                pbestFit(rand10(i), :)       = n_fit(pareto(j), :);
                break;
            end
        end
        % gbest
        t_p1_chrom = [gbest_p1_chrom; n_p1_chrom(pareto, :)];
        t_m1_chrom = [gbest_m1_chrom; n_m1_chrom(pareto, :)];
        t_p2_chrom = [gbest_p2_chrom; n_p2_chrom(pareto, :)];
        t_m2_chrom = [gbest_m2_chrom; n_m2_chrom(pareto, :)];
        t_fit      = [gbestFit; n_fit(pareto, :)];
        [t_fit, idx] = unique(t_fit, 'rows'); % ШЅжи
        t_p1_chrom = t_p1_chrom(idx, :);
        t_m1_chrom = t_m1_chrom(idx, :);
        t_p2_chrom = t_p2_chrom(idx, :);
        t_m2_chrom = t_m2_chrom(idx, :);
        [~, pareto] = nonRank(t_fit);
        gbest_p1_chrom = t_p1_chrom(pareto, :);
        gbest_m1_chrom = t_m1_chrom(pareto, :);
        gbest_p2_chrom = t_p2_chrom(pareto, :);
        gbest_m2_chrom = t_m2_chrom(pareto, :);
        gbestFit       = t_fit(pareto, :);
        % tbest
        [Min, idx] = min(n_fit(:, 1));
        if Min < tbestFit(1, 1)
            tbest_p1_chrom = n_p1_chrom(idx, :);
            tbest_m1_chrom = n_m1_chrom(idx, :);
            tbest_p2_chrom = n_p2_chrom(idx, :);
            tbest_m2_chrom = n_m2_chrom(idx, :);
            tbestFit       = n_fit(idx, :);
        end
        % cbest
        [Min, idx] = min(n_fit(:, 2));
        if Min < cbestFit(1, 2)
            cbest_p1_chrom = n_p1_chrom(idx, :);
            cbest_m1_chrom = n_m1_chrom(idx, :);
            cbest_p2_chrom = n_p2_chrom(idx, :);
            cbest_m2_chrom = n_m2_chrom(idx, :);
            cbestFit       = n_fit(idx, :);
        end
    end
    %% local search tbestFit cbestFit
    localSearchBest(1);
    localSearchBest(2);
    idx = randperm(PS, 2);
    p1_chrom(idx(1), :) = tbest_p1_chrom;
    m1_chrom(idx(1), :) = tbest_m1_chrom;
    p2_chrom(idx(1), :) = tbest_p2_chrom;
    m2_chrom(idx(1), :) = tbest_m2_chrom;
    fit(idx(1), :)      = tbestFit;
    p1_chrom(idx(2), :) = cbest_p1_chrom;
    m1_chrom(idx(2), :) = cbest_m1_chrom;
    p2_chrom(idx(2), :) = cbest_p2_chrom;
    m2_chrom(idx(2), :) = cbest_m2_chrom;
    fit(idx(2), :)      = cbestFit;
    %% anti-stagnation
    if (iter > 0.1 * Itermax && size(gbestFit, 1) == preNum)
        initNum = ceil(PS/10);
        [new_p1_chrom, new_m1_chrom, new_p2_chrom, new_m2_chrom] = initPop(initNum);
        new_fit = zeros(initNum, 2);
        for i = 1 : initNum
            new_fit(i, :) = calFit(new_p1_chrom(i, :), new_m1_chrom(i, :), new_p2_chrom(i, :), new_m2_chrom(i, :));
        end
        idx = randperm(PS, initNum);
        p1_chrom(idx, :) = new_p1_chrom;
        m1_chrom(idx, :) = new_m1_chrom;
        p2_chrom(idx, :) = new_p2_chrom;
        m2_chrom(idx, :) = new_m2_chrom;
        fit(idx, :)      = new_fit;
    end
    preNum = size(gbestFit, 1);
end
%% elite
elite = [gbest_p1_chrom gbest_p2_chrom gbest_m1_chrom gbest_m2_chrom];
elite_fit = gbestFit;