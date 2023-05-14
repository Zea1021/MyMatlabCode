function indiv = insert(indiv, pos)
pos = sort(pos);
temp1 = indiv(1 : pos(2) - 1);
temp2 = indiv(pos(2) : end);
temp = temp1(pos(1));
temp1(pos(1))=[];
indiv=[temp1,temp,temp2];