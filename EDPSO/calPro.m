function [finish, totalEner, mt] = calPro(p2_chrom, m2_chrom, finish_d, totalEner)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3;

mt = zeros(1, M2);      % completion time of machine
tm = cell(M2, 1);       % machine time record
finish_p = [];          % completion time of repro

s1 = p2_chrom;  % operation genes of repro
p = ones(N, C); % the next operation

for i = 1 : SH
    signal = s1(i);
    job = floor(signal/100); 
    compo = signal - job * 100; 
    k = p(job, compo); % operation
    p(job, compo) = p(job, compo) + 1; % update the next
    m = m2_chrom(i); % machine
    
    if (k == 1)
        avaiTime = finish_d(job);
    else
        avaiTime = finish_p(job, compo, k - 1);
    end
    
    if (~isempty(tm{m})) % find gap
        Ts = zeros(1, size(tm{m}, 2)); % init gap start
        Te = zeros(1, size(tm{m}, 2)); % init gap end
        Ts(1) = avaiTime; 
        Te(1) = tm{m}(1, 1); 
        for it = 1 : size(tm{m}, 2) - 1
            Ts(it + 1) = max([avaiTime, tm{m}(2, it)]);
            Te(it + 1) = tm{m}(1, it + 1);
        end
        gasp = Te - Ts;
        pos = find(gasp >= Time2(job, compo, k, m), 1); % judge gap
        if ~isempty(pos) % meeting
            startTime = Ts(pos);
            endTime = startTime + Time2(job, compo, k, m);
            % update mt{m}
            tm{m} = [tm{m}(:, 1 : pos-1), [startTime; endTime], tm{m}(:, pos : end)]; % ¸üÐÂtm
            % disp('left-shift success');
        else             % not meeting
            if (mt(m) < avaiTime) % init startTime
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
    finish_p(job, compo, k) = endTime;
    totalEner = totalEner + CM2(m) * Time2(job, compo, k, m);
end

finish = [];
for i = 1 : N
    for j = 1 : C
        finish(i, j) = finish_p(i, j, H(i, j)); % the repro completion time of Cij
    end
end