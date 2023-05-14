function SubChrom=aberrance1(SubChrom,Pm,F)
global n p k MD MP ct Time total_L total_p;
num_pop=size(SubChrom,1);
for i=1:num_pop
    if rand<Pm           
        val=SubChrom(i,:);
        if rand<F
            pos=randi(n,1,1);             
            val=aberrance_minT(val,pos);
        else
            pos=randi(total_L-n,1,1)+n;   
            val=aberrance_minT(val,pos);
        end
        SubChrom(i,:)=val;
    end
end
