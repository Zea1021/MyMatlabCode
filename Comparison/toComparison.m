function [BestMeanObj, meanC_metric, MeanVariIGD, MeanVariHV] = toComparison(CaseNo)
clearvars -except CaseNo;
%% pre
prefix = 'data_';
suffix = '_result_';
AllInstance = {'4-4', '4-5', '6-4', '6-5', '8-8', '8-10', '10-8', '10-10'; ...
    'P4-C16-O72', 'P4-C20-O76', 'P6-C24-O84', 'P6-C30-O108', ...
    'P8-C64-O216', 'P8-C80-O296', 'P10-C80-O230', 'P10-C100-O360'};
data = AllInstance{1, CaseNo};
num = 20;
file = {'IMOGA/', 'NSGA-II/', 'PSOGA/', 'MDFFO/', 'MODE/', 'mMOEAD/', 'EDPSO/'};
algorithm = {'IMOGA', 'NSGA-II', 'PSOGA', 'MDFFO', 'MODE', 'm-MOEA/D', 'EDPSO'};
algoNum = length(algorithm);
%% init
truePF = [];
minVal = [inf inf];
maxVal = [0 0];
%% load data
objBest = cell(algoNum, 1);
for k = 1 : algoNum
    for i = 1 : num
        load([file{k} prefix data suffix num2str(i)]);
        truePF(end + 1 : end + size(elite_fit, 1), :) = elite_fit;
        truePF = unique(truePF, 'rows'); 
        [~, pareto] = nonRank(truePF);
        truePF = truePF(pareto, :);
        if size(elite_fit, 1) == 1
            minVal = min(minVal, elite_fit);
            maxVal = max(maxVal, elite_fit);
            objBest{k}(i, :) = elite_fit;
        else
            minVal = min(minVal, min(elite_fit));
            maxVal = max(maxVal, max(elite_fit));
            objBest{k}(i, :) = min(elite_fit);
        end
    end
end
BestMeanObj = zeros(2, 2 * algoNum);
for k = 1 : algoNum
    for i = 1 : 2
        BestMeanObj(i, 1 + 2 * (k - 1)) = min(objBest{k}(:, i));
        BestMeanObj(i, 2 + 2 * (k - 1)) = mean(objBest{k}(:, i));
    end
end
%% C_metric
C_metric = cell(1, algoNum - 1);
meanC_metric = zeros(1, 2 * (algoNum - 1));
for k = 2 : algoNum
    for i = 1 : num
        load([file{1} prefix data suffix num2str(i)]);
        temp = elite_fit;
        load([file{k} prefix data suffix num2str(i)]);
        C_metric{k - 1}(i, :) = calC_metric(temp, elite_fit);
    end
    meanC_metric(1, 1 + 2 * (k - 2) : 2 + 2 * (k - 2)) = mean(C_metric{k - 1});
end
%% IGD HV
IGD = zeros(num, algoNum);
HV = zeros(num, algoNum);
for k = 1 : algoNum
    for i = 1 : num
        load([file{k} prefix data suffix num2str(i)]);
        IGD(i, k) = calIGD(elite_fit, truePF, minVal, maxVal);
        HV(i, k) = calHV(elite_fit, minVal, maxVal);
    end
end
MeanVariIGD = zeros(1, 2 * algoNum);
MeanVariHV = zeros(1, 2 * algoNum);
for k = 1 : algoNum
    MeanVariIGD(1, 1 + 2 * (k - 1)) = mean(IGD(:, k));
    MeanVariIGD(1, 2 + 2 * (k - 1)) = std(IGD(:, k));
    MeanVariHV(1, 1 + 2 * (k - 1))  = mean(HV(:, k));
    MeanVariHV(1, 2 + 2 * (k - 1))  = std(HV(:, k));
end
%% best idx of IGD
bestIdx = zeros(algoNum, 1);
for k = 1 : algoNum
    [~, idx] = min(IGD(:, k));
    bestIdx(k) = idx;
end