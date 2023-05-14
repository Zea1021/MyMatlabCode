function chrom = initMethod1(original, Pi)
global n p k MD MP ct Time total_L total_p;
m1_chrom=zeros(1,n);
m2_chrom=zeros(1,total_p);
%% dis
p1_chrom = randperm(n);
for i = 1 : n
    job = p1_chrom(1, i);
    if rand <= Pi
        [~, m] = min(Time(job).Dis);
    else
        m = randi(MD);
    end
    m1_chrom(1, i) = m;
end
%% repro
p2_chrom = original(randperm(total_p));
nextPro = ones(n, p); % the next operation
for i = 1 : total_p
    signal = p2_chrom(i);
    job = floor(signal/100);   % product
    comp = signal - job * 100; % component
    pro = nextPro(job, comp);  % operation
    nextPro(job, comp) = nextPro(job, comp) + 1; % update
    if rand <= Pi
        [~, m] = min(Time(job).Pro{comp}(pro, :));
    else
        m = randi(MP);
    end
    m2_chrom(1, i) = m;
end
chrom = [p1_chrom p2_chrom m1_chrom m2_chrom];