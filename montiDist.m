function p = montiDist(AllC,Iter)

% function p = montiDist(AllC,Iter)
% This function proposes a  probability distribution of simulaneous activity
% by randomizing the entries in each column of AllC and summing, Iter
% times. Does not preserve autostructure of individual responses.

% it is a good idea to seed rand if that has not been done this session:
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
% seed from clock.

% Finn Upham 2012-04-20
 
AltC = AllC;
L = size(AllC);

v = zeros(1,L(2)+1);

for k = 1:Iter
    [~,RandI] = sort(rand(size(AllC)));

    for i = 1:L(2)
        AltC(:,i)=AllC(RandI(:,i),i);
    end
    
    v = v + hist(sum(AltC,2),0:L(2))/L(1);
    
end

p = v'/Iter;
