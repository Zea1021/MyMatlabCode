function Pa = calImproDegree(Mean)
Pa = zeros(1, size(Mean, 2) - 1);
for i = 2 : size(Mean, 2)
    Pa(i - 1) = abs(Mean(1) - Mean(i))/(10 *  Mean(i));
end
