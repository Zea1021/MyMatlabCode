function Pop = initPopByRule(PS,ruleNum)
global n p k MD MP ct Time total_L total_p;
Pop=zeros(PS, 2 * total_L); 
original = zeros(1, total_p);
count = 1;
for i = 1 : n
    for j = 1 : p
        for l = 1 : k(j)
            original(count) = i * 100 + j;
            count = count + 1;
        end
    end
end
for i = 1 : ruleNum
    Pop(i, :) = initMethod1(original, 0.5);
end
for i = ruleNum + 1 : 2 * ruleNum
    Pop(i, :) = initMethod2(original, 0.5);
end
for i = 2 * ruleNum + 1 : PS
    p1_chrom = randperm(n);
    p2_chrom = original(randperm(total_p));
    m1_chrom = randi(MD,1,n);
    m2_chrom = randi(MP,1,total_p);
    Pop(i, :) = [p1_chrom p2_chrom m1_chrom m2_chrom];
end