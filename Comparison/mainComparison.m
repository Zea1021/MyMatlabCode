clear
close all
clc

Obj = [];
C_metric = [];
IGD = [];
HV = [];
for i = 1 : 8
    [BestMeanObj, meanC_metric, MeanVariIGD, MeanVariHV] = toComparison(i);
    Obj(end+1:end+2, :)      = BestMeanObj;
    C_metric(end+1:end+1, :) = meanC_metric;
    IGD(end+1:end+1, :)      = MeanVariIGD;
    HV(end+1:end+1, :)       = MeanVariHV;
end
mIGD = mean(IGD);
mHV = mean(HV);