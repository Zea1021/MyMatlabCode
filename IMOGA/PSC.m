% PSC
function newChrom = PSC(chrom)
global n p k MD MP ct Time total_L total_p;
p1_chrom = chrom(1 : n); 
p2_chrom = chrom(n+1 : total_L);
m1_chrom = chrom(total_L+1 : total_L+n);
m2_chrom = chrom(total_L+n+1 : 2 * total_L);
new_p1_chrom = p1_chrom;
new_m1_chrom = m1_chrom;
new_p2_chrom = p2_chrom;
new_m2_chrom = m2_chrom;

p_gene = chrom(1 : total_L);
m_gene = chrom(total_L + 1 : 2 * total_L);
%% dis
S = p_gene(1, 1 : n);  % operation 
M = m_gene(1, 1 : n);  % machine
finish_d = zeros(2,n); % completion time of each product
mt = zeros(1, MD);     % completion time of each machine
for i = 1 : n
    job = S(i);
    m = M(i);  
    startTime = mt(m);
    endTime = startTime + Time(job).Dis(m);
    mt(m) = endTime;
    finish_d(job) = endTime;
end
%% repeo
S = p_gene(n + 1 : end);  
M = m_gene(n + 1 : end);  
finish_c = zeros(n, p, max(k)); % completion time of Cijk
mt = zeros(1, MP);       % completion time of machine
tm = cell(MP, 1);        % machine time record
nextPro = ones(n, p);    % the next operation of Cij
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
    
    flag = false; % not left shift
    if (~isempty(tm{m})) % find gap
        Ts = zeros(1, size(tm{m}, 2)); % init
        Te = zeros(1, size(tm{m}, 2)); 
        Ts(1) = avaiTime; 
        Te(1) = tm{m}(1, 1); 
        for it = 1 : size(tm{m}, 2) - 1
            Ts(it + 1) = max([avaiTime, tm{m}(2, it)]);
            Te(it + 1) = tm{m}(1, it + 1);
        end
        gasp = Te - Ts;
        pos = find(gasp >= Time(job).Pro{comp}(pro, m), 1); % judge gap
        if ~isempty(pos) % meeting
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
end
finish_p = zeros(n, p);
for i = 1 : n
    for j = 1 : p
        finish_p(i, j) = finish_c(i, j, k(j)); 
    end
end
%% Ass
finish = finish_p;
record = zeros(n, p);
for i = 1 : n
    [~, compJ] = min(finish);
    record(i, :) = compJ;
    for j = 1 : p 
        finish(compJ(j), j) = inf; 
    end
end
%% change
job = randi(n); % remanufactured product
compJ = record(job, :); 
set = zeros(1, p);  % component set
temp = zeros(1, p); % time
for i = 1 : p
    temp(i) = finish_p(compJ(i), i);
    set(i) = compJ(i) * 100 + i;
end
%% adjust
[~, idx] = max(temp); 
signal = compJ(idx) * 100 + idx;
idx = find(p2_chrom == signal, 1 ,'last'); 

temp = randi(total_p);
while (ismember(p2_chrom(temp), set))
    temp = randi(total_p);
end
new_p2_chrom([idx temp]) = new_p2_chrom([temp idx]);
new_m2_chrom([idx temp]) = new_m2_chrom([temp idx]);
newChrom = [new_p1_chrom new_p2_chrom new_m1_chrom new_m2_chrom];