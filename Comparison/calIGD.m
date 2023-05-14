function IGD = calIGD(PF, truePF, minVals, maxVals)
m1 = size(PF, 1);
m = size(truePF, 1);
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