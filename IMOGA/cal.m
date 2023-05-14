function fit = cal(chrom)
global n p k MD MP ct Time total_L;

p_gene = chrom(1 : total_L);
m_gene = chrom(total_L + 1 : 2 * total_L);
energy = 0;
%% dis
S = p_gene(1, 1 : n);  % operation
M = m_gene(1, 1 : n);  % machine
finish_d = zeros(2,n); % init dis completion time
mt = zeros(1, MD);     % init machine time
for i = 1 : n
    job = S(i);% product
    m = M(i);  % machine
    startTime = mt(m);
    endTime = startTime + Time(job).Dis(m);
    mt(m) = endTime;
    energy = energy + ct{1}(m) * Time(job).Dis(m);
    finish_d(job) = endTime;
end
%% repro
S = p_gene(n + 1 : end);  
M = m_gene(n + 1 : end); 
finish_c = zeros(n, p, max(k)); % completion time of Cijk
mt = zeros(1, MP);       % init machine time
tm = cell(MP, 1);        % machine time record
nextPro = ones(n, p);    % the next repro operation of Cij
for i = 1 : length(S)
    signal = S(i);
    job = floor(signal/100);   % product
    comp = signal - job * 100; % component
    pro = nextPro(job, comp);  % operation of Cij
    nextPro(job, comp) = nextPro(job, comp) + 1; % update
    m = M(i);                   
    
    if (pro == 1)
        avaiTime = finish_d(job);
    else
        avaiTime = finish_c(job, comp, pro - 1);
    end
    
    flag = false; % not left shift
    if (~isempty(tm{m})) % find gap
        Ts = zeros(1, size(tm{m}, 2)); % init gap startTime and endTime
        Te = zeros(1, size(tm{m}, 2)); 
        Ts(1) = avaiTime; 
        Te(1) = tm{m}(1, 1); 
        for it = 1 : size(tm{m}, 2) - 1
            Ts(it + 1) = max([avaiTime, tm{m}(2, it)]);
            Te(it + 1) = tm{m}(1, it + 1);
        end
        gasp = Te - Ts;
        pos = find(gasp >= Time(job).Pro{comp}(pro, m), 1); % judge gap
        if ~isempty(pos) 
            startTime = Ts(pos);
            endTime = startTime + Time(job).Pro{comp}(pro, m);
            tm{m} = [tm{m}(:, 1 : pos-1), [startTime; endTime], tm{m}(:, pos : end)]; % update tm
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
finish_p = zeros(n, p);
for i = 1 : n
    for j = 1 : p
        finish_p(i, j) = finish_c(i, j, k(j)); % the repro completion time of Cij
    end
end
%% Ass
finish = finish_p;
m = 1;
mt = zeros(1, m); % init machine time
record = zeros(n, p);
for i = 1 : n
    [minTime, compJ] = min(finish);
    record(i, :) = compJ;
    Max = max(minTime);
    if(mt(m) < Max)
        startTime = Max;
    else
        startTime = mt(m); 
    end
    endTime = startTime + Time(i).Ass; 
    mt(m) = endTime; 
    energy = energy + ct{3}(m) * Time(i).Ass; 
    
    for j = 1 : p 
        finish(compJ(j), j) = inf; 
    end
end
finish_a = max(mt);
fit = [finish_a energy];