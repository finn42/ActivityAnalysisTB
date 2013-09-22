function p = empDist(Coinc,empLike)

x = [0:max(max([Coinc empLike]))];

rH = hist(Coinc,x)';%/length(Coinc(Time>0,:));

L = size(empLike);

V = zeros(length(x),L(2));

for i = 1:L(2)
    V(:,i) = hist(empLike(:,i),x)'; %/L(1);
end

Vmean = mean(V,2);

dMA = sum((rH-Vmean).^2)^0.5;

dMV = sum((V-repmat(Vmean,1,L(2))).^2,1).^0.5;

p = length(dMV(dMV>dMA))/L(2);