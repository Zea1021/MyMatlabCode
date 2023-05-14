function val1=aberrance_minT(val1,pos)
global n p k MD MP ct Time total_L total_p;
if pos<=n
    num_n=val1(1,pos);             
    T=Time(num_n).Dis;
    [~,index]=min(T);
    val1(1,total_L+pos)=index;     
else 
    num_n=floor(val1(pos)/100);
    num_p=val1(1,pos)-num_n*100;
    index=find(val1==val1(pos));
    num_k= index==pos;
    T=Time(num_n).Pro{num_p}(num_k,:);
    [~,index]=min(T);
    val1(1,total_L+pos)=index;
end

