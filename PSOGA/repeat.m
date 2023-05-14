clear;
clc;

repeatNum = 20;
Data = {'4-4','4-5','6-4','6-5','8-8','8-10','10-8','10-10'};
prefix = 'data_';
suffix = '_result_';
saveRoad = 'result/';
if ~exist(saveRoad,'dir')
    mkdir(saveRoad);
end
for i = 1 : size(Data, 2)
    filename = Data{i};
    for j = 1 : repeatNum
        [elite, elite_fit] = mainLoop(filename);
        save([saveRoad, [prefix filename suffix num2str(j)]], 'elite', 'elite_fit');
    end
end