% input:
%        PF:自己的前沿面
%        truePF：真实前沿面，由所有的解求出
%        minVals:目标函数的最小值
%        maxVals:目标函数的最大值
% output:
%        IGD:IGD评价指标

% [~, Sort] = sort(PF, 'ascend');
% PF = PF(Sort, :);
% plot(PF(:, 1), PF(:, 2),'o-b')
% hold on;
% [~, Sort] = sort(truePF, 'ascend');
% truePF = truePF(Sort, :);
% plot(truePF(:, 1), truePF(:, 2),'*-r')
% hold on;
function IGD = calIGD(PF, truePF, minVals, maxVals)
m1 = size(PF, 1);
m = size(truePF, 1);
% Get the normalized front
PF = (PF - repmat(minVals, m1, 1)) ./ repmat(maxVals - minVals, m1, 1);
truePF = (truePF - repmat(minVals, m, 1)) ./ repmat(maxVals - minVals, m, 1);
IGD = 0;
for i = 1:m
    diff = repmat(truePF(i,:), m1, 1) - PF;
    dist = sqrt(sum(diff.^2, 2));
    IGD = IGD + min(dist);
end
IGD = IGD/m;

end


% function IGD = IGD_matlab(PF, truePF)
% 
% q = 2; %the parameter of IGD
% %STEP 1. Obtain the maximum and minimum values of the Pareto front
% m1 = size(PF, 1);
% m = size(truePF, 1);
% maxVals = max(truePF);
% minVals = min(truePF);
% 
% %STEP 2. Get the normalized front
% normalizedPF = (PF - repmat(minVals, m1, 1)) ./ repmat(maxVals - minVals, m1, 1);
% normalizedTruePF = (truePF - repmat(minVals, m, 1)) ./ repmat(maxVals - minVals, m, 1);
% 
% %STEP 3. Sum the distances between each point of the front and the nearest point in the true Pareto front
% IGD = 0;
% for i = 1:m
%     diff = repmat(normalizedTruePF(i,:), m1, 1) - normalizedPF;
%     dist = sqrt(sum(diff.^2, 2));
%     
%     IGD = IGD + min(dist)^q;
% end
% IGD = IGD^(1.0/q)/m;
% 
% end