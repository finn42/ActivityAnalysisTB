function [cSa,cSm] = monoCoord(Time,Series,winSize,Thresh,option,nonPara)

Iter = 200;

cSa = zeros(winSize,1);
cSm = cSa;

N = size(Series,2);
binN = 4;

if nargin < 6
    nonPara = 1;
end

for i=1:winSize

    [AC,AllC]=actionCount(Time(i:end),Series(i:end,:),winSize,winSize,Thresh,option);
    [~,cSa(i)]=sichiSq(AC,N,binN);
    
    if nonPara == 1
        p = montiDist(AllC,Iter);
        [~,cSm(i)]=sichiSq(AC,N,binN,p);
    end
end
 
cSa=mean(-log10(cSa+10^(-16)));
cSm=mean(-log10(cSm+10^(-16)));