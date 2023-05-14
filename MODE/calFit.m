function Fit = calFit(Pop)
Size = size(Pop, 1);
Fit = zeros(Size, 2);
for i = 1 : Size
   Fit(i, :) = cal(Pop(i, :));
end