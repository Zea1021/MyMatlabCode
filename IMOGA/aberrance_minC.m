function val=aberrance_minC(val,pos)
global n p k MD MP ct Time total_L total_p;
if pos<=n
    job = val(pos);
    temp = ct{1}(1) * Time(job).Dis(1);
    idx = 1;
    for i = 2 : MD
        if ct{1}(i) * Time(job).Dis(i) < temp
            temp = ct{1}(i) * Time(job).Dis(i);
            idx = i;
        end
    end
    val(1,total_L+pos) = idx;
else
    signal = val(pos);
    job = floor(signal/100); 
    comp = signal - job * 100; 
    idx = find(val == signal);
    pro = find(idx == pos); 
    temp = ct{2}(1) * Time(job).Pro{comp}(pro, 1);
    idx = 1;
    for i = 2 : MP
        if ct{2}(i) * Time(job).Pro{comp}(pro, i) < temp
            temp = ct{2}(i) * Time(job).Pro{comp}(pro, i);
            idx = i;
        end
    end
    val(1,total_L+pos) = idx;
end