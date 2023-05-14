function indiv = levorotation(indiv, pos)
pos = sort(pos);
temp=indiv(pos(2));
indiv(pos(2))=[];
indiv=[indiv(1:pos(1)-1),temp,indiv(pos(1):end)];