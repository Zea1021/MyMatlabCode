function Box(Algorithm, Data)
algoNum = length(Algorithm);
X = [];
Y = {};
for algo = 1 : algoNum
    X(end + 1 : end + size(Data, 1)) = Data(:, algo);
    Y(end + 1 : end + size(Data, 1)) = repmat(Algorithm(algo), size(Data, 1), 1);
end
boxplot(X, Y);