function [C,pval,DAct,bins]=sichiSq(AC,N,k,P)

%function [C,pval,DAct,bins]=sichiSq(AC,N,k,P);

% this function performes a goodness of fit on the activity distribution AC
% with an optimisation process to selecting appropriate bins according to
% the binomial random model generated as the null hypothesis. Note AC is
% only the activity count time series, a single column without time
% stamps. N is the number of participants in the population that generated
% AC and k is the number of bins to use in the calculation. 

% If a fourth argument, P, is included as the probability distrubution for
% activity levels of 0:N/N, in a column.

% The output Dact is a two column matrix with the second being the actual
% activity histogram for all numbers of participants (0:N), and the second
% being the estimated distribution of the binomial model.
% bins is two columns as well (model, actual) of the bins selected for the
% goodness of fit test.
% C is the chi squared value for the two set of bins. 

% pval is the likelyhood of a distribution as extreme or more so being
% generated by the binomial model, each as summarized in the bins.
% This function requires the Statistics Toolbox and equiSplit

% Finn Upham, March 31st, 2011
% updated 2012/08/22

L=size(AC,1);


% to calculate the null hypothesis of a binomial random model
% avg rate of activity

n=N;
aL = (0:n)';
p = zeros(N+1,1);

if nargin==4
    % use the subplied activity probability to set up the random model.
    if length(P)==1
        a = P;
        for i=0:n
            p(i+1)=((a^i)*(1-a)^(n-i))*...
                (factorial(n)/(factorial(n-i)*factorial(i)));
        end
    elseif length(P)==N+1
        p = P;
        
        if size(p,2)~=1
            p = p';
        end
    else
       error('MatLab:sichiSq:DistModel',...
           'Distribution must be of suitable form, 1x1 or (N+1)x1') 
    end
    
    dof = 0;
    
else
    a=sum(AC)/(L);
    % calculate the binomial probability of different activity levels.
    for i=0:n
        p(i+1)=((a^i)*(1-a)^(n-i))*...
            (factorial(n)/(factorial(n-i)*factorial(i)));
    end
    dof = -1;
end


pAct=L*p;
NAct=hist(AC*N,aL);
DAct=[NAct' pAct];

[v,bins] = equiSplit(pAct,k,5);

k = size(bins,1);

bins = [zeros(size(bins)) bins];

for i = 1:k
    bins(i,1) = sum(NAct(v{i}));
end

C=sum(sum(((bins(:,1)-bins(:,2)).^2)./bins(:,1)));

pval=1-chi2cdf(C,k-1+dof);

% if isnan(pval)
%     pval = 1;
% end






 