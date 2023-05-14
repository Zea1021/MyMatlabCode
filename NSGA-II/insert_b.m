function indiv=insert_b(indiv,position)
position=sort(position);
temp=indiv(position(2));
indiv(position(2))=[];
indiv=[indiv(1:position(1)-1),temp,indiv(position(1):end)];