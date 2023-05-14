function C_metric = calC_metric(PF1,PF2)
C_metric=zeros(1,2);
num1=size(PF1,1);
num2=size(PF2,1);

num=0;
for i=1:num2
    val2=PF2(i,:);
    flag=0;
    for j=1:num1
        val1=PF1(j,:);
        if (val1(1,1)<=val2(1,1)&&val1(1,2)<val2(1,2))||(val1(1,1)<val2(1,1)&&val1(1,2)<=val2(1,2))
            flag=1;
            break;
        end
    end
    if flag==1
        num=num+1;
    end
end
C_metric(1,1)=num/num2;

num=0;
for i=1:num1
    val1=PF1(i,:);
    flag=0;
    for j=1:num2
        val2=PF2(j,:);
        if (val2(1,1)<=val1(1,1)&&val2(1,2)<val1(1,2))||(val2(1,1)<val1(1,1)&&val2(1,2)<=val1(1,2))
            flag=1;
            break;
        end
    end
    if flag==1
        num=num+1;
    end
end
C_metric(1,2)=num/num1;