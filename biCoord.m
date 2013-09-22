function [cSa] = biCoord(Time,SeriesA,SeriesB,winSize,Thresh,option1,option2)

cSa = zeros(winSize,1);

if nargin < 7
    option2 = option1;
end

for i=1:winSize

    [~,AllCa]=actionCount(Time(i:end),SeriesA(i:end,:),winSize,...
        winSize,Thresh,option1);
    [~,AllCb]=actionCount(Time(i:end),SeriesB(i:end,:),winSize,...
        winSize,Thresh,option2);
    
    [~,cSa(i)]=jointChiSq(AllCa,AllCb,'Ind');
end

cSa=mean(-log10(cSa+10^(-16)));