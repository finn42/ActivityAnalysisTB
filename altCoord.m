function [cSa] = altCoord(Time,Series,winSize,ThreshA,optionA,ThreshB,optionB)

cSa = zeros(winSize,1);

N = size(Series,2);

for i=1:winSize

    [~,AllCa]=actionCount(Time(i:end),Series(i:end,:),winSize,...
        winSize,ThreshA,optionA);
    [~,AllCb]=actionCount(Time(i:end),Series(i:end,:),winSize,...
        winSize,ThreshB,optionB);
    
    [~,cSa(i)]=jointChiSq(AllCa,AllCb,'Alt');
end

cSa=mean(-log10(cSa+10^(-16)));