function Subval=toCross(val1,val2)
global n p k MD MP ct Time total_L total_p;
%% operation
val1Chunk1=val1(1:n);
val2Chunk1=val2(1:n);
val1Chunk2=val1(n+1:total_L);
val2Chunk2=val2(n+1:total_L);
[val1Chunk1sub,val2Chunk1sub]=poxcross(val1Chunk1,val2Chunk1);
[val1Chunk2sub,val2Chunk2sub]=poxcross(val1Chunk2,val2Chunk2);
%% machine
val1SM=val1(total_L+1:end);
val2SM=val2(total_L+1:end);
val1SM1=val1SM(1:n);
val2SM1=val2SM(1:n);
val1SM2=val1SM(n+1:total_L);
val2SM2=val2SM(n+1:total_L);
[val1SM1sub,val2SM1sub]=duodiancross(val1SM1,val2SM1);
[val1SM2sub,val2SM2sub]=duodiancross(val1SM2,val2SM2);
%% record
Subval=zeros(2,2*total_L);
Subval(1,:)=[val1Chunk1sub,val1Chunk2sub,val1SM1sub,val1SM2sub];
Subval(2,:)=[val2Chunk1sub,val2Chunk2sub,val2SM1sub,val2SM2sub];