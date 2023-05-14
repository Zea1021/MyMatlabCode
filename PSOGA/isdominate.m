function flag = isdominate(fit1, fit2)
dom_less=0;
dom_equal=0;
dom_more=0;
m = size(fit1, 2);
for i = 1 : m % for each obj
    if (fit1(i) < fit2(i))
        dom_less = dom_less + 1;
    elseif (fit1(i) == fit2(i))
        dom_equal = dom_equal + 1;
    else
        dom_more = dom_more + 1;
    end
end
if dom_less == 0 && dom_equal ~= m % fit1 is dominated by fit2 
    flag = -1;
elseif dom_more == 0 && dom_equal ~= m % fit1 dominates fit2 
    flag = 1;
else
    flag = 0;
end