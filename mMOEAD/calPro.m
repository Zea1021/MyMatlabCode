function [finish_p, totalEner, mt, finish_c] = calPro(p2_chrom, m2_chrom, finish_d, totalEner)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3;

mt = zeros(1, M2);      
tm = cell(M2, 1);          
finish_c = zeros(N, C, max(H(:))); 

s1 = p2_chrom;
p = ones(N, C); 

for i = 1 : SH
    signal = s1(i);
    job = floor(signal/100); 
    compo = signal - job * 100; 
    k = p(job, compo); 
    p(job, compo) = p(job, compo) + 1; 
    m = m2_chrom(i); 
    
    if (k == 1)
        avaiTime = finish_d(job); 
    else
        avaiTime = finish_c(job, compo, k - 1);
    end
    
    flag = false; 
    if (~isempty(tm{m})) 
        Ts = zeros(1, size(tm{m}, 2)); 
        Te = zeros(1, size(tm{m}, 2)); 
        Ts(1) = avaiTime; 
        Te(1) = tm{m}(1, 1); 
        for it = 1 : size(tm{m}, 2) - 1
            Ts(it + 1) = max([avaiTime, tm{m}(2, it)]);
            Te(it + 1) = tm{m}(1, it + 1);
        end
        gasp = Te - Ts;
        pos = find(gasp >= Time2(job, compo, k, m), 1);
        if ~isempty(pos) 
            flag = true;
            startTime = Ts(pos);
            endTime = startTime + Time2(job, compo, k, m);
            tm{m} = [tm{m}(:, 1 : pos-1), [startTime; endTime], tm{m}(:, pos : end)]; 
        end
    end
    if (~flag) 
        if (mt(m) < avaiTime) 
            startTime = avaiTime; 
        else 
            startTime = mt(m);
        end
        endTime = startTime + Time2(job, compo, k, m);
        mt(m) = endTime;
        tm{m} = [tm{m} [startTime; endTime]];
    end
    finish_c(job, compo, k) = endTime;
    totalEner = totalEner + CM2(m) * Time2(job, compo, k, m); 
end

finish_p = zeros(N, C);
for i = 1 : N
    for j = 1 : C
        finish_p(i, j) = finish_c(i, j, H(i, j));
    end
end