function [n_p1_chrom, n_m1_chrom, n_p2_chrom, n_m2_chrom, n_fit] = localSearch(p1_chrom, m1_chrom, p2_chrom, m2_chrom, fit)
localTime = 10;
count = 0;
n_p1_chrom = zeros(localTime, size(p1_chrom, 2));
n_m1_chrom = zeros(localTime, size(p1_chrom, 2));
n_p2_chrom = zeros(localTime, size(p2_chrom, 2));
n_m2_chrom = zeros(localTime, size(p2_chrom, 2));
n_fit      = zeros(localTime, 2);
while count < localTime
    % levorotation
    pos = randperm(length(p1_chrom), 2);
    t_p1_chrom = levorotation(p1_chrom, pos);
    t_m1_chrom = levorotation(m1_chrom, pos);
    pos = randperm(length(p2_chrom), 2);
    t_p2_chrom = levorotation(p2_chrom, pos);
    t_m2_chrom = levorotation(m2_chrom, pos);
    % cal fit
    t_fit = calFit(t_p1_chrom, t_m1_chrom, t_p2_chrom, t_m2_chrom);
    flag = isdominate(t_fit, fit);
    if flag ~= 1
        pos = randperm(length(t_p1_chrom), 2);
        t_p1_chrom = insert(t_p1_chrom, pos);
        t_m1_chrom = insert(t_m1_chrom, pos);
        pos = randperm(length(t_p2_chrom), 2);
        t_p2_chrom = insert(t_p2_chrom, pos);
        t_m2_chrom = insert(t_m2_chrom, pos);
        t_fit = calFit(t_p1_chrom, t_m1_chrom, t_p2_chrom, t_m2_chrom);
    end
    
    n_p1_chrom(count + 1, :) = t_p1_chrom;
    n_m1_chrom(count + 1, :) = t_m1_chrom;
    n_p2_chrom(count + 1, :) = t_p2_chrom;
    n_m2_chrom(count + 1, :) = t_m2_chrom;
    n_fit(count + 1, :) = t_fit;
    count = count + 1;
end