function [AC,AllC,dT]=actionCount(Time,Series,FrameSize,HopSize,Thresh,option)

%function [AC,AllC,dT]=actionCount(Time,Series,FrameSize,HopSize,Thresh,option)
%
% Version 3.0 (earlier versions are compatible with older scripts
%
% function to assess the rate of activity for the columns of "Series" over
% frames of size FrameSize + 1 sample at intervals of Hopsize. 
% Action option must be one of:
% 'Change','Inc', 'Dec','UBound','Percent','LMax','LMin','UBound','Xup' 
% If the type of activity specified by option is found in a column to a
% degree exceeding the threshold 'thresh', that column is marked as active
% (1) for that time frame centered at the time cited in corresponding row 
% of the dT series. Not LMax and LMin not optimal at present.

% AC reports the proportion of colums showing the specified activity in
% each time frame.

% Finn Upham 2012 09 06
 
if nargin==5
    option='Change';
end

if size(Time,1)~=size(Series,1)
error('Matlab:actionCount:RowsMatch',...
        'Rows of Time and Series must match')
end

% set up the frames for activity analysis and define the time points
% associated with each frame.

k=HopSize;
l=FrameSize;

if k>l
error('Matlab:actionCount2:HopSize',...
        'Warning hopsize larger than framesize, may miss activity.')
end

T=Time(round(k/2):k:size(Series,1),:);
dT=T;
frame=cell(length(T),1);
cSize = size(Series,2);

for i=1:length(T)
    if k*(i-1)+l>size(Series,1)
        frame{i}=Series(k*(i-1):end,:);
    elseif i ==1
        frame{i}=Series(1:l,:);
    else
        frame{i}=Series(k*(i-1):k*(i-1)+l,:);
    end
end

%Now evaluate and save activity over each frame, according to the requested
%activity type

AC=zeros(size(T));
AllC=zeros(size(T,1),size(Series,2));

if strcmp(option, 'Inc')
    for i=1:length(T)
        %take the difference of the first and last frame    
        D=frame{i}(end,:)-frame{i}(1,:);
        % Set D(r)==1 if response r increased over frame i at least Thresh,
        % else D(r)==0
        D(D<Thresh)=0;
        D=sign(D);
        %Save the individual activity and activity count for frame i
        AllC(i,:)=D;
        AC(i)=sum(D)/cSize;
    end
    
elseif strcmp(option, 'Dec')
    for i=1:length(T)
        %take the difference of the first and last frame    
        D=frame{i}(end,:)-frame{i}(1,:);
        %Set D(r)==1 if response r decreased over frame i, else D(r)==0
        D(D>-Thresh)=0;
        D=abs(sign(D));
        %Save the individual activity and activity count for frame i
        AllC(i,:)=D;
        AC(i)=sum(D)/cSize;
    end
    
elseif strcmp(option, 'Change')
    for i=1:length(T)
        %take the difference of the first and last frame    
        D=abs(frame{i}(end,:)-frame{i}(1,:));
        %Set D(r)==1 if response r changes over frame i, else D(r)==0
        D(D<Thresh)=0;
        D=abs(sign(D));
        %Save the individual activity and activity count for frame i
        AllC(i,:)=D;
        AC(i)=sum(D)/cSize;
    end
elseif strcmp(option, 'UBound')
    for i=1:length(T)
        %Set the UBound to zero
        D=frame{i};
        D=D-Thresh;
        %Push all non-negative response values to 1, and rest to 0
        D(D>=0)=1;
        D(D<0)=0;
        %sum over frame to catch which response have any points with values
        %over or at the UBound and set their activity values to 1 for the frame
        D=sign(sum(D,1));
        %Save the individual activity and activity count for frame i
        AllC(i,:)=D;
        AC(i)=sum(D)/cSize;
    end
elseif strcmp(option, 'LBound')
    for i=1:length(T)
        %Set the LBound to zero
        D=frame{i};
        D=D-Thresh;
        %Push all non-positive response values to 1, and rest to 0
        D(D<=0)=0;
        D(D>0)=1;
        D = 1-D;
        %sum over frame to catch which response have any points with values
        %below or at the LBound and set their activity values to 1 for the frame
        D=sign(sum(D,1));
        %Save the individual activity and activity count for frame i
        AllC(i,:)=D;
        AC(i)=sum(D)/cSize;
    end
elseif strcmp(option, 'Percent')
    %to identify the frames with the highest average values (in place of
    %sample wise rank value), we first average response values over each
    %frame
    
    V=AllC;
    for i=1:length(T)
        V(i,:)=mean(frame{i});
    end
    if Thresh<0
        V=-V;
        Thresh=-Thresh;
    end
    %find the percentil values for each response over these frame averages
    Y=prctile(V,Thresh);
    
    %throwing out instances when the percentile includes too much of the
    %signal
    for i=1:size(V,2)
        if length(V(V(:,i)>=Y(i)))>length(T)*2*(100-Thresh)/100
            V(:,i)=Y(i)-1;
        end
    end
    
    for i=1:length(T)
        %per frame, treat the percentiles values as UBounds per responses
        D=V(i,:)-Y;
        %Push all non-negative response values to 1, and rest to 0
        D(D>=0)=1;
        D(D<0)=0;
        %Save the individual activity and activity count for frame i
        AllC(i,:)=D;
        AC(i)=sum(D)/cSize;
    end
elseif strcmp(option, 'LMax')
    for i=1:length(T)
        %Set the LBound to zero
        if k<l
            D=frame{i};
        elseif k==l
            if i < length(T)
                D=[frame{i};frame{i+1}(1:2,:)];
            end
        end
        D2 = D;
        D=diff([D(1,:); D; D(end,:)],1,1);
        D = diff(sign(D),1,1);
        D(D>0)=0;
        D(D2<Thresh)=0;
        D=abs(sign(sum(D,1)));
        %Save the individual activity and activity count for frame i
        AllC(i,:)=D;
        AC(i)=sum(D)/cSize;
    end
elseif strcmp(option, 'LMin')
    for i=1:length(T)
        %Set the LBound to zero
        if k<l
            D=frame{i};
        elseif k==l
            if i < length(T)
                D=[frame{i};frame{i+1}(1,:)];
            end
        end
        D=diff(D,1,1);
        D = diff(sign(D),1,1);
        D(D<0)=0;
        D=abs(sign(sum(D,1)));
        %Save the individual activity and activity count for frame i
        AllC(i,:)=D;
        AC(i)=sum(D)/cSize;
    end
elseif strcmp(option, 'Xup')
    for i=1:length(T)
        if k<l
            D=frame{i};
        elseif k==l
            if i < length(T)
                D=[frame{i};frame{i+1}(1,:)];
            end
        end
        %Set the threshold to zero to capture crossing upwards
        D = sign(D-Thresh);
        D=diff(D,1,1);
        D(D<0)=0;
        D=abs(sign(sum(D,1)));
        %Save the individual activity and activity count for frame i
        AllC(i,:)=D;
        AC(i)=sum(D)/cSize;
    end
else
     error('MATLAB:actionCount:ppOutput', ...
            'option not supported')

end




