clear;
close all;
clc;

global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3; 
global PS Pc Pm;
%% data
load DataSet/10-10.mat
PS = 120;  
Pc = 0.9; 
Pm = 0.2; 
T  = 10;  
maxIter = 500; 
%% init weight and neighbour
[W,B] = Init_weights(PS - 1,2, T);
W(W==0) = 0.000001;
%% init pop
ruleNum = ceil(PS/10);
[p1_chrom, p2_chrom, m1_chrom, m2_chrom, fit] = initPopByRule(PS,ruleNum);
%% refer point
Z = min(fit);
%% main loop
for iter = 1 : maxIter
    for i = 1 : PS
        if rand < Pc
            P = B(i, :); % neighbour
            k = randperm(length(P));
            [c_p1_chrom, c_m1_chrom, c_p2_chrom, c_m2_chrom] = evolution(p1_chrom(P(k(1)),:), m1_chrom(P(k(1)),:), p2_chrom(P(k(1)),:), m2_chrom(P(k(1)),:), p1_chrom(P(k(2)),:), m1_chrom(P(k(2)),:), p2_chrom(P(k(2)),:), m2_chrom(P(k(2)),:));
            c_fit = calFit(c_p1_chrom, c_m1_chrom, c_p2_chrom, c_m2_chrom);
            
            % update Z
            Z = min(Z, c_fit);
            % update neighbour
            for j = randperm(length(P))
                g_old = max(abs(fit(P(j), :) - Z) .* W(P(j), :));
                g_new = max(abs(c_fit - Z) .* W(P(j), :));
                if g_new < g_old
                    p1_chrom(P(j),:) = c_p1_chrom;
                    m1_chrom(P(j),:) = c_m1_chrom;
                    p2_chrom(P(j),:) = c_p2_chrom;
                    m2_chrom(P(j),:) = c_m2_chrom;
                    fit(P(j),:)      = c_fit;
                end
            end
        end
    end
end
%% elite
[unique_fit, idx] = unique(fit, 'rows'); 
unique_p1_chrom = p1_chrom(idx, :);
unique_m1_chrom = m1_chrom(idx, :);
unique_p2_chrom = p2_chrom(idx, :);
unique_m2_chrom = m2_chrom(idx, :);
[~, pareto]     = nonRank(unique_fit);
elite_p1_chrom = unique_p1_chrom(pareto, :);
elite_m1_chrom = unique_m1_chrom(pareto, :);
elite_p2_chrom = unique_p2_chrom(pareto, :);
elite_m2_chrom = unique_m2_chrom(pareto, :);
elite = [elite_p1_chrom elite_p2_chrom elite_m1_chrom elite_m2_chrom];
elite_fit = unique_fit(pareto, :);
