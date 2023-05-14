function chrom = initMethod2(original, Pi)
global n p k MD MP ct Time total_L total_p;
m1_chrom=zeros(1,n);
m2_chrom=zeros(1,total_p);
%% diss
p1_chrom = randperm(n);
for i = 1 : n
    if rand <= Pi
        [~, m] = min(ct{1});
    else
        m = randi(MD);
    end
    m1_chrom(1, i) = m;
end
%% repro
p2_chrom = original(randperm(total_p));
nextPro = ones(n, p); 
for i = 1 : total_p
    signal = p2_chrom(i);
    job = floor(signal/100); 
    comp = signal - job * 100; 
    pro = nextPro(job, comp); 
    nextPro(job, comp) = nextPro(job, comp) + 1; 
    if rand <= Pi
        m = 1;
        for j = 2 : MP
            if (ct{2}(j) * Time(job).Pro{comp}(pro, j) < ct{2}(m) * Time(job).Pro{comp}(pro, m))
                m = j;
            end
        end
    else
        m = randi(MP);
    end
    m2_chrom(1, i) = m;
end
chrom = [p1_chrom p2_chrom m1_chrom m2_chrom];