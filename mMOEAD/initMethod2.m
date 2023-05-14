function [p1_chrom, m1_chrom, p2_chrom, m2_chrom, fit] = initMethod2(Pi)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3;
m1_chrom=zeros(1,N);
m2_chrom=zeros(1,sum(sum(H)));

original = zeros(1, SH);
count = 1;
for j = 1 : N
    for k = 1 : C
        for l = 1 : H(j, k)
            original(count) = j * 100 + k;
            count = count + 1;
        end
    end
end

%% dis
p1_chrom = randperm(N);

mt = zeros(1, M1);      
finish_d = zeros(1, N); 
totalEner = 0; 
for i = 1 : N
    job = p1_chrom(1, i);
    if rand <= Pi
        m = 1;
        for j = 2 : M1
            if (Time1(job, j) < Time1(job, m))
                m = j;
            end
        end
    else
        m = randi(M1, 1, 1);
    end
    m1_chrom(1, i) = m;
    
    startTime = mt(m);
    endTime = startTime + Time1(job, m);
    mt(m) = endTime;  
    finish_d(job) = endTime;           
    totalEner = totalEner + CM1(m) * Time1(job, m); 
end

%% repro
p2_chrom = original(randperm(SH));

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
    if rand <= Pi
        m = 1;
        for j = 2 : M2
            if (Time2(job, compo, k, j) < Time2(job, compo, k, m))
                m = j;
            end
        end
    else
        m = randi(M2, 1, 1);
    end
    m2_chrom(1, i) = m;
    
    if (k == 1)
        avaiTime = finish_d(job); 
    else
        avaiTime = finish_c(job, compo, k - 1);
    end
    
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
            startTime = Ts(pos);
            endTime = startTime + Time2(job, compo, k, m);
            tm{m} = [tm{m}(:, 1 : pos-1), [startTime; endTime], tm{m}(:, pos : end)];
        else           
            if (mt(m) < avaiTime) 
                startTime = avaiTime; 
            else 
                startTime = mt(m);
            end
            endTime = startTime + Time2(job, compo, k, m);
            mt(m) = endTime;
            tm{m} = [tm{m} [startTime; endTime]];
        end
    else
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

%% ass
[finish_a, totalEner] = calAss(finish_p, totalEner);

fit = [finish_a totalEner];