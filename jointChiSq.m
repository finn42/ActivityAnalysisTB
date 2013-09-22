function [C,p,DAct,Bins,v1,v2]=jointChiSq(AllC1,AllC2,option,k)
% function [C,p,DAct,Bins,v1,v2]=jointChiSq(AllC1,AllC2,option,k)
% This function performes a goodness of fit test on the joint activity 
% distribution of two activities, assessed on responses over matching time 
% frames against the assumption that these are independent. AllC* are
% matricies, one column per responses in each collection, each row the time
% frame, entries 1 if response j is active in time frame i, else 0. 

% also enter option 'Ind' for activitys which are independe, say from 
% same activity from different groups of responses or forms of activity 
% which can happen simultaneously in the same response
% or 'Alt' for alterating activity, like rating increases and decreases,
% which necessarily alternate by exploring the other. The random model is
% generated differently for these two circumstances.

% The output Dact is a 3D array with the first matrix being the actual
% joint distribution for all numbers of participants, and the second
% being the estimated distribution if the series were independent.
% Bins is 3X3X2 array as well (actual, model) of the bins selected for the
% goodness of fit test. If the bins hold too few samples (< 5), the bins 
% are recalculated to generate 6 bins rather than 9 to describe the
% distribution.

% C is the chi squared value from the two sets of bins. 

% p is the likelyhood of the actual distribution being the result of
% independent processes. This requires the Statistics Toolbox. Also depends
% on function equiSlit and DistAlt.

% Finn Upham, April 6th, 2012
% Updated 2012/08/23 Alternating solved! 

C = 0;
p = 0;
DAct = [];
Bins = [];
v1 = [];
v2 = [];

AC1 = sum(AllC1,2);
N1 = size(AllC1,2);

AC2 = sum(AllC2,2);
N2 = size(AllC2,2);

if nargin<4
    k = 3;
end

L=size(AC1,1);

if size(AC2,1)~=L
    error('MatLab:jointChiSq',...
        'The activity time series are not of matching length')
end



% Calcualte the joint distributions, actual and model.

% independent distribution, for independent activity
 
if strcmp(option,'Ind')
    % actual distribution
    A=zeros(N1+1,N2+1);
    for i=0:N1
        A(i+1,:)=hist(AC2(AC1==i),0:N2); %make sure this isn't backwards
    end
    
    % calculate the three bins for distribution of AC1
    [NAct1]=hist(AC1,0:N1);
    pref = 'Even';
    [v1,bins1] = equiSplit(NAct1,k,16,pref);

    % calculate the three bins for distribution of AC2
    [NAct2]=hist(AC2,0:N2);
    [v2,bins2] = equiSplit(NAct2,k,16,pref);
    
    % expected distribution over all joint activity-levels
    B = NAct1'*NAct2/L;
    
    % dof = 0; % degrees of freedom adjustment
    
    % actual distribution
    A=zeros(N1+1,N2+1);
    for j=0:N1
        A(j+1,:)=hist(AC2(AC1==j),0:N2); %make sure this isn't backwards
    end
end

% alternating distribution
if strcmp(option,'Alt')
    % check if actually alternating
%     if (sum(AllC1(logical(AllC2)))+sum(AllC2(logical(AllC1))))>0
%         error('Activities 1 and 2 must be exclusive to apply this test.')
%     end
    
    % derive hypothetical activity from total activity
    [A,B] = JointDistAlt(AllC1,AllC2);
    pref = 'Even';
    NAct1 = sum(B,2);
    NAct2 = sum(B,1);
    
    [v1,bins1] = equiSplit(NAct1,k,16,pref);
    [v2,bins2] = equiSplit(NAct2,k,16,pref);
    % dof = -1; % one (or two) estimated variables in alt distribution
    
end


% In case the binnings of each separate dimension did not yeild sufficient
% expected samples in each table entry, we reduce the number of rows or
% columns accordingly.

%check out the models table to see if the binnings in both dimensions
%resulted in appropriate distribution of samples in all entries of the
%table
tableB = zeros(length(bins1),length(bins2));

if isempty(tableB)
    C = 0;
    p = 1;
    DAct = [];
    Bins = cell(0);
    v1 = [];
    v2 = [];
    return;
end
    
for i = 1:length(bins1)
    for j = 1:length(bins2)
        tableB(i,j) = sum(sum(B(v1{i},v2{j})));
    end
end

%check the quality of the model distribution, shrink the table if bad.
while min(min(tableB))<5
    
    if max(size(tableB))>2
        
        if length(bins2)<length(bins1)
            [v1,bins1] = equiSplit(NAct1,k-1,16,pref);
        elseif length(bins2)>length(bins1)
            [v2,bins2] = equiSplit(NAct2,k-1,16,pref);
        elseif  min(bins1)<min(bins2)
            [v1,bins1] = equiSplit(NAct1,k-1,16,pref);
            k = k-1;
        else
            [v2,bins2] = equiSplit(NAct2,k-1,16,pref);
            k = k-1;
        end

        tableB = zeros(length(bins1),length(bins2));
        for i = 1:length(bins1)
            for j = 1:length(bins2)
                tableB(i,j) = sum(sum(B(v1{i},v2{j})));
            end
        end
    else
        fprintf('Bad distributions, can not test these activity levels/n')
        p = NaN;
        C = NaN;
        return;
    end
end


% table for actual
tableA = zeros(size(tableB));
for i = 1:length(bins1)
    for j = 1:length(bins2)
        tableA(i,j) = sum(sum(A(v1{i},v2{j})));
    end
end


DAct(:,:,1)=B;
DAct(:,:,2)=A;

% Finally calculate the Goodness of fit test

Bins=zeros(length(bins1),length(bins2),2);
Bins(:,:,1)=tableA;
Bins(:,:,2)=tableB;

C=sum(sum(((Bins(:,:,2)-Bins(:,:,1)).^2)./Bins(:,:,2)));

p=1-chi2cdf(C,length(bins1)+length(bins2)-2);



 