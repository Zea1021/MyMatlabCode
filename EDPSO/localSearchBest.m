function localSearchBest(code)
global tbest_p1_chrom tbest_m1_chrom tbest_p2_chrom tbest_m2_chrom tbestFit;
global cbest_p1_chrom cbest_m1_chrom cbest_p2_chrom cbest_m2_chrom cbestFit;
localTime = 10;
count = 0;
while count < localTime
    if code == 1
        p1_chrom = tbest_p1_chrom;
        m1_chrom = tbest_m1_chrom;
        p2_chrom = tbest_p2_chrom;
        m2_chrom = tbest_m2_chrom;
        fit      = tbestFit;
    else
        p1_chrom = cbest_p1_chrom;
        m1_chrom = cbest_m1_chrom;
        p2_chrom = cbest_p2_chrom;
        m2_chrom = cbest_m2_chrom;
        fit      = cbestFit;
    end
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
    count = count + 1;
    
    % tbest
    if t_fit(1, 1) < tbestFit(1, 1)
        tbest_p1_chrom = t_p1_chrom;
        tbest_m1_chrom = t_m1_chrom;
        tbest_p2_chrom = t_p2_chrom;
        tbest_m2_chrom = t_m2_chrom;
        tbestFit = t_fit;
    end
    % cbest
    if t_fit(1, 2) < cbestFit(1, 2)
        cbest_p1_chrom = t_p1_chrom;
        cbest_m1_chrom = t_m1_chrom;
        cbest_p2_chrom = t_p2_chrom;
        cbest_m2_chrom = t_m2_chrom;
        cbestFit = t_fit;
    end
end

