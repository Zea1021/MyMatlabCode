function newChrom = BWL(chrom, F)
global n p k MD MP ct Time total_L;

p1_chrom = chrom(1 : n); 
p2_chrom = chrom(n+1 : total_L);
m1_chrom = chrom(total_L+1 : total_L+n);
m2_chrom = chrom(total_L+n+1 : 2 * total_L);
new_p1_chrom = p1_chrom;
new_m1_chrom = m1_chrom;
new_p2_chrom = p2_chrom;
new_m2_chrom = m2_chrom;
if rand < F
    result = histc(m1_chrom, 1 : MD);
    minVal = min(result);                         
    lowM = find(result == minVal);  
    maxVal = max(result);
    highM = find(result == maxVal); 
    
    mt = zeros(1, MD);      
    for i = 1 : n
        job = p1_chrom(1, i);
        m = m1_chrom(1, i);
        mt(m) = mt(m) + Time(job).Dis(m);
    end
    
    [~, minIdx] = min(mt);  
    [~, maxIdx] = max(mt);  
    if (ismember(minIdx, lowM) && ismember(maxIdx, highM))
        pos = find(m1_chrom == maxIdx);
        pos = pos(randi(length(pos)));
        new_m1_chrom(pos) = minIdx;
    else
        highM = highM(randi(length(highM)));
        result(result == 0) = 0.01;
        pp = cumsum(1 ./ (result)) / sum(1 ./ (result));
        idx = find(rand <= pp, 1);
        while (idx == highM)
            idx = find(rand <= pp, 1);
        end
        pos = find(m1_chrom == highM);
        pos = pos(randi(length(pos)));
        new_m1_chrom(pos) = idx;
    end
else
    result = histc(m2_chrom, 1 : MP);
    minVal = min(result);                         
    lowM = find(result == minVal); 
    maxVal = max(result);
    highM = find(result == maxVal);
    
    p_gene = chrom(1 : total_L);
    m_gene = chrom(total_L + 1 : 2 * total_L);
    energy = 0;
    %% dis
    S = p_gene(1, 1 : n);  
    M = m_gene(1, 1 : n);  
    finish_d = zeros(2,n); 
    mt = zeros(1, MD);  
    for i = 1 : n
        job = S(i);
        m = M(i);  
        startTime = mt(m);
        endTime = startTime + Time(job).Dis(m);
        mt(m) = endTime;
        energy = energy + ct{1}(m) * Time(job).Dis(m);
        finish_d(job) = endTime;
    end
    %% repro
    S = p_gene(n + 1 : end);  
    M = m_gene(n + 1 : end);  
    finish_c = zeros(n, p, max(k)); 
    mt = zeros(1, MP);       
    tm = cell(MP, 1);       
    nextPro = ones(n, p);    
    for i = 1 : length(S)
        signal = S(i);
        job = floor(signal/100);    
        comp = signal - job * 100; 
        pro = nextPro(job, comp); 
        nextPro(job, comp) = nextPro(job, comp) + 1; 
        m = M(i);                   
        
        if (pro == 1)
            avaiTime = finish_d(job); 
        else
            avaiTime = finish_c(job, comp, pro - 1);
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
            pos = find(gasp >= Time(job).Pro{comp}(pro, m), 1); 
            if ~isempty(pos) 
                startTime = Ts(pos);
                endTime = startTime + Time(job).Pro{comp}(pro, m);
                tm{m} = [tm{m}(:, 1 : pos-1), [startTime; endTime], tm{m}(:, pos : end)]; 
                flag = true; 
            end
        end
        if (~flag) 
            if (mt(m) < avaiTime) 
                startTime = avaiTime; 
            else 
                startTime = mt(m);
            end
            endTime = startTime + Time(job).Pro{comp}(pro, m);
            mt(m) = endTime;
            tm{m} = [tm{m} [startTime; endTime]];
        end
        finish_c(job, comp, pro) = endTime;
        energy = energy + ct{2}(m) * Time(job).Pro{comp}(pro, m);
    end
    
    [~, minIdx] = min(mt);  
    [~, maxIdx] = max(mt);  
    if (ismember(minIdx, lowM) && ismember(maxIdx, highM))
        pos = find(m2_chrom == maxIdx);
        pos = pos(randi(length(pos)));
        new_m2_chrom(pos) = minIdx;
    else
        highM = highM(randi(length(highM)));
        result(result == 0) = 0.01;
        pp = cumsum(1 ./ (result)) / sum(1 ./ (result));
        idx = find(rand <= pp, 1);
        while (idx == highM)
            idx = find(rand <= pp, 1);
        end
        pos = find(m2_chrom == highM);
        pos = pos(randi(length(pos)));
        new_m2_chrom(pos) = idx;
    end
end
newChrom = [new_p1_chrom new_p2_chrom new_m1_chrom new_m2_chrom];