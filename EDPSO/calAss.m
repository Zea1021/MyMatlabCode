function [finish_a, totalEner, record] = calAss(finish_p, totalEner)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3;

mt = zeros(1, M3); % completion time of machine
finish = finish_p;
record = zeros(N, C); % component selection
for i = 1 : N
    [minTime, idx1] = min(finish); 
    record(i, :) = idx1;
    [~, idx2] = max(minTime);
    job = idx1(idx2);
    compo = idx2;
    m = 1;
    if(mt(m) < finish_p(job, compo)) 
        mt(m) = finish_p(job, compo) + Time3(m); 
    else 
        mt(m) = mt(m) + Time3(m);
    end
    totalEner = totalEner + CM3(m) * Time3(m); 
    
    for j = 1 : C
        finish(idx1(j), j) = inf;
    end
end

finish_a = mt(1);
for i = 2 : M3
    if(mt(i) < finish_a)
        finish_a = mt(i);
    end
end






