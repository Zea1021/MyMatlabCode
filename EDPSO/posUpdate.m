function [c_p1_chrom, c_m1_chrom, c_p2_chrom, c_m2_chrom, c_fit] = posUpdate(p1_chrom, m1_chrom, p2_chrom, m2_chrom)
global N C H SH M1 M2 M3 CM1 CM2 CM3 Time1 Time2 Time3; 
global Pp Pg Pt
global pbest_p1_chrom pbest_m1_chrom pbest_p2_chrom pbest_m2_chrom pbestFit;
global gbest_p1_chrom gbest_m1_chrom gbest_p2_chrom gbest_m2_chrom gbestFit;
global tbest_p1_chrom tbest_m1_chrom tbest_p2_chrom tbest_m2_chrom tbestFit;
global cbest_p1_chrom cbest_m1_chrom cbest_p2_chrom cbest_m2_chrom cbestFit;

Size = size(p1_chrom, 1);
c_p1_chrom = [];
c_m1_chrom = [];
c_p2_chrom = [];
c_m2_chrom = [];
c_fit = [];
for i = 1 : Size
    t1_p1_chrom = p1_chrom(i, :);
    t1_m1_chrom = m1_chrom(i, :);
    t1_p2_chrom = p2_chrom(i, :);
    t1_m2_chrom = m2_chrom(i, :);
    %% swap
    pos = randperm(N, 2);
    t1_p1_chrom([pos(2), pos(1)]) = t1_p1_chrom([pos(1), pos(2)]);
    t1_m1_chrom([pos(2), pos(1)]) = t1_m1_chrom([pos(1), pos(2)]);
    t1_p2_chrom([pos(2), pos(1)]) = t1_p2_chrom([pos(1), pos(2)]);
    t1_m2_chrom([pos(2), pos(1)]) = t1_m2_chrom([pos(1), pos(2)]);
    
    %% cross
    r = rand;
    if r < Pp
        t2_p1_chrom = pbest_p1_chrom(i, :);
        t2_m1_chrom = pbest_m1_chrom(i, :);
        t2_p2_chrom = pbest_p2_chrom(i, :);
        t2_m2_chrom = pbest_m2_chrom(i, :);
    elseif r < Pp + Pg
        idx = randi(size(gbest_p1_chrom, 1));
        t2_p1_chrom = gbest_p1_chrom(idx, :);
        t2_m1_chrom = gbest_m1_chrom(idx, :);
        t2_p2_chrom = gbest_p2_chrom(idx, :);
        t2_m2_chrom = gbest_m2_chrom(idx, :);
    elseif r < Pp + Pg + Pt
        t2_p1_chrom = tbest_p1_chrom;
        t2_m1_chrom = tbest_m1_chrom;
        t2_p2_chrom = tbest_p2_chrom;
        t2_m2_chrom = tbest_m2_chrom;
    else
        t2_p1_chrom = cbest_p1_chrom;
        t2_m1_chrom = cbest_m1_chrom;
        t2_p2_chrom = cbest_p2_chrom;
        t2_m2_chrom = cbest_m2_chrom;
    end
    
    if rand < 0.5
        [new_p1_chrom1, new_p1_chrom2] = poxCross(t1_p1_chrom, t2_p1_chrom);
        [new_p2_chrom1, new_p2_chrom2] = poxCross(t1_p2_chrom, t2_p2_chrom);
    else
        [new_p1_chrom1, new_p2_chrom1, new_p1_chrom2, new_p2_chrom2] = preCross(t1_p1_chrom, t1_p2_chrom, t2_p1_chrom, t2_p2_chrom);
    end
    [new_m1_chrom1, new_m1_chrom2] = binarycross(t1_m1_chrom, t2_m1_chrom);
    [new_m2_chrom1, new_m2_chrom2] = binarycross(t1_m2_chrom, t2_m2_chrom);
    
    % cal fit
    new_fit1 = calFit(new_p1_chrom1, new_m1_chrom1, new_p2_chrom1, new_m2_chrom1);
    new_fit2 = calFit(new_p1_chrom2, new_m1_chrom2, new_p2_chrom2, new_m2_chrom2);
    
    % pbest
    flag = isdominate(new_fit1, pbestFit(i, :));
    if flag == 1
        pbest_p1_chrom(i, :) = new_p1_chrom1;
        pbest_m1_chrom(i, :) = new_m1_chrom1;
        pbest_p2_chrom(i, :) = new_p2_chrom1;
        pbest_m2_chrom(i, :) = new_m2_chrom1;
        pbestFit(i, :) = new_fit1;
    end
    flag = isdominate(new_fit2, pbestFit(i, :));
    if flag == 1
        pbest_p1_chrom(i, :) = new_p1_chrom2;
        pbest_m1_chrom(i, :) = new_m1_chrom2;
        pbest_p2_chrom(i, :) = new_p2_chrom2;
        pbest_m2_chrom(i, :) = new_m2_chrom2;
        pbestFit(i, :) = new_fit2;
    end
    
    c_p1_chrom(end + 1 : end + 2, :) = [new_p1_chrom1; new_p1_chrom2];
    c_m1_chrom(end + 1 : end + 2, :) = [new_m1_chrom1; new_m1_chrom2];
    c_p2_chrom(end + 1 : end + 2, :) = [new_p2_chrom1; new_p2_chrom2];
    c_m2_chrom(end + 1 : end + 2, :) = [new_m2_chrom1; new_m2_chrom2];
    c_fit(end + 1 : end + 2, :)      = [new_fit1; new_fit2];
end
% tbest
if new_fit2(1, 1) < tbestFit(1, 1)
    tbest_p1_chrom = new_p1_chrom2;
    tbest_m1_chrom = new_m1_chrom2;
    tbest_p2_chrom = new_p2_chrom2;
    tbest_m2_chrom = new_m2_chrom2;
    tbestFit = new_fit2;
end
% cbest
if new_fit2(1, 2) < cbestFit(1, 2)
    cbest_p1_chrom = new_p1_chrom2;
    cbest_m1_chrom = new_m1_chrom2;
    cbest_p2_chrom = new_p2_chrom2;
    cbest_m2_chrom = new_m2_chrom2;
    cbestFit = new_fit2;
end