function HV = calHV(PF,minVals,maxVals)
numPoint  = size(PF,1);
PF=sortrows([PF(:,1) PF(:,2)],[1 -2]);
PF = (PF - repmat(minVals, numPoint, 1)) ./ repmat(maxVals - minVals, numPoint, 1);
X=zeros(1,numPoint);
Y=zeros(1,numPoint);
for i=1:numPoint
    if i==numPoint
        X(i)=1-PF(i,1);
        Y(i)=1-PF(i,2);
        break;
    end
    X(i)=PF(i+1,1)-PF(i,1);
    Y(i)=1-PF(i,2);
end
HV=sum(X.*Y);