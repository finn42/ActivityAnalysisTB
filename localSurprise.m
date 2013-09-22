function [Coinc,ECoinc,pCoinc,sCoinc,empLike,sigCoinc]=localSurprise(data,tw,tr,N,option)

% function [Coinc,ECoinc,pCoinc,sCoinc,empLike,sigCoinc]=LocalSurprise(data,tw,tr,N,option)
% Version 1.0
% localSurprise assess the likelyhood of the coincidence of events in a
% collection of time synchronised point processes in windows of size tw,
% by non-parametrically estimating the distribution of coincidences by
% repeatedly randomly shifting the alignment of the series uniformly over
% time span tr. 
% Inputs
% data: a collection of point process time series {0,1}, one series per column,
%       sampled at an constant sample rate and mutually aligned 
% tw: the window of coincidence, in samples
% tr: the window of random shifting, in samples
% N: the interations of random shifting to inform the expecte distribution
% option: either 'Loop' or 'Extend' this specifies the method of
%   approximating the expecte distribution for the ends of the series. Loop,
%   the default, pulls data from the opposite end of the series rather than
%   zero pad. Extend only evaluates the ranges which are not effected by the
%   ends points and then eextends the first and alst distributions to conver
%   each extremity.
%
% Outputs
% Coinc: a time series reporting the number of time series active in each
%       time frame, essentially the sum of data across columns.
% ECoinc: the mean of the estimated distribution for each time frame
% pCoinc: the p value or percentile of the actual Coinc count per time
%   frame in the estimated distribution.
% sCoinc: the p-values translated into surprise values (-log10((1-p)/p))
% emplike: the coincidences per random shuffle in each time frame.
% sigCoinc: the local coincidence threshold for p<0.01, per time frame
%
% All of these outputs have the same number of rows as data.

% this method of estimating local surprise is derived from:
% Pipa, G., Wheeler, D., Singer, W., and Nikoli ?c, D. (2008). Neuroxidence:
%   reliable and efficient analysis of an excess or deficiency of joint-
%   spike events. Journal of computational neuroscience, 25(1):64?88.

% Finn Upham 2013 07 16

if nargin==4
    option='Loop';
end

L = size(data);

% building alt data
empLike = zeros(L(1),N);

for j = 1:N
    shifts = round((rand(1,L(2))-0.5)*tr);
    
    Sdata = zeros(L+[tr 0]);

    for i = 1:L(2)
        d = data(:,i);
        % pad the begining and end of the shifted data with the
        % opposite end of the series
        Sdata(:,i) = [d(end-((tr/2)+shifts(i)-1):end,1); data(:,i);d(1:((tr/2)-shifts(i)),1)];
    end
    Sdata = Sdata((tr/2)+1:end-(tr/2),:);
    empLike(:,j) = sum(sign(conv2(ones(tw,1),1,Sdata,'same')),2); 
end

n =  sum(sign(conv2(ones(tw,1),1,data,'same')),2);

p = zeros(size(n));
q = p;

if strcmp(option,'Extend')
    empLike = [repmat(empLike(1+tr/2,:),tr/2,1);...
        empLike((tr/2)+1:end-(tr/2),:);repmat(empLike(end-tr/2,:),tr/2,1)];
end

for i = 1:length(p)
    dist = empLike(i,:);
    p(i) = max([length(dist(dist>n(i))),1])/N;
    q(i) = prctile(dist,99);
end

popSur = log10((1-p)./p);
popSur(popSur>3) = 3;
popSur(popSur<-3) = -3;


Coinc = n;
ECoinc = mean(empLike,2);
pCoinc = p;
sCoinc = popSur;
sigCoinc = q;
