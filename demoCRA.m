%% load a data set
clear

load Behavioural.mat
% Coll{1} = Beh{19}.Data;
% Time{1} = Beh{19}.Time;
% Tag{1} = [Beh{19}.Piece ' ' Beh{19}.Audience ' ' Beh{19}.Measure];
% 
% Coll{2} = Beh{20}.Data;
% Time{2} = Beh{20}.Time;
% Tag{2} = [Beh{20}.Piece ' ' Beh{20}.Audience ' ' Beh{20}.Measure];
% 
% Coll{3} = Beh{23}.Data;
% Time{3} = Beh{23}.Time;
% Tag{3} = [Beh{23}.Piece ' ' Beh{23}.Audience ' ' Beh{23}.Measure];
% 
% 
% Coll{4} = Beh{24}.Data;
% Time{4} = Beh{24}.Time;
% Tag{4} = [Beh{24}.Piece ' ' Beh{24}.Audience ' ' Beh{24}.Measure];

Coll{1} = Beh{21}.Data;
Time{1} = Beh{21}.Time;
Tag{1} = [Beh{21}.Piece ' ' Beh{21}.Audience ' ' Beh{21}.Measure];

Coll{2} = Beh{22}.Data;
Time{2} = Beh{22}.Time;
Tag{2} = [Beh{22}.Piece ' ' Beh{22}.Audience ' ' Beh{22}.Measure];

Coll{3} = Beh{25}.Data;
Time{3} = Beh{25}.Time;
Tag{3} = [Beh{25}.Piece ' ' Beh{25}.Audience ' ' Beh{25}.Measure];


Coll{4} = Beh{26}.Data;
Time{4} = Beh{26}.Time;
Tag{4} = [Beh{26}.Piece ' ' Beh{26}.Audience ' ' Beh{26}.Measure];

%% plot the basics of activity
i = 1;
Data = Coll{i}(:,1:3);
T = Time{i};

figure
subplot(4,1,1)
plot(T,Data)
axis([T(1) T(end) 0 1])
    xlabel('Time (s)')
    ylabel('Rating Scale')
    title([Tag{i} ' a few responses'])
    
subplot(4,1,2)
[ACinc,AllC,dT] = actionCount(T,Data,4,4,0.05,'Inc');
stem(dT,AllC)
axis tight
xlabel('Time (s)')
 ylabel('Increasing activity')
 
 subplot(4,1,3)
[ACinc,AllC,dT] = actionCount(T,Data,4,4,0.05,'Dec');
stem(dT,AllC)
axis tight
xlabel('Time (s)')
 ylabel('Decreasing activity')
 
  subplot(4,1,4)
[ACinc,AllC,dT] = actionCount(T,Data,4,4,0.05,'LMin');
stem(dT,AllC)
axis tight
xlabel('Time (s)')
 ylabel('Local Minima activity')

%% plot the data

figure('WindowStyle','docked')

for i = 1:4
    subplot(4,2,(i-1)*2 +1)
    plot(Time{i},Coll{i})
    axis([Time{i}(1) Time{i}(end) 0 1])
    xlabel('Time (s)')
    ylabel('Rating Scale')
    title([Tag{i} ' all responses'])

    subplot(4,2,(i-1)*2 +2)
    [ACinc,~,dT] = actionCount(Time{i},Coll{i},4,4,0.05,'Inc');
    ACdec = actionCount(Time{i},Coll{i},4,4,0.05,'Dec');
    stem(dT,[ACinc -ACdec])
    axis([Time{i}(1) Time{i}(end) -1 1])
    legend Increases Decreases
    xlabel('Time (s)')
    ylabel('Activity-level')
    title([Tag{i} ' rating change activity'])
end

%% test the activities first simple

% increasing activity
figure('WindowStyle','docked')
for i = 1:4
    subplot(4,3,[1 2]+(i-1)*3)
    [ACinc,AllC,dT] = actionCount(Time{i},Coll{i},4,4,0.05,'Inc');
    stem(dT,ACinc)
    axis([Time{i}(1) Time{i}(end) 0 1])
    xlabel('Time (s)')
    ylabel('Activity-level')
    title([Tag{i} ' rating increases activity'])
    
    [C,pval,DAct,bins]=sichiSq(ACinc,size(AllC,2),5);
    subplot(4,6,5+(i-1)*6)
    bar([0:size(AllC,2)]/size(AllC,2),DAct)
    legend('Collection','Model','Location','NorthEastOutside')
    ylabel('# timeframes')
    xlabel('Activity-level')
    axis tight
    title('Activity-level distributions')
    
    subplot(4,6,6+(i-1)*6)
    bar(bins)
    legend('Collection','Model','Location','NorthEastOutside')
    ylabel('# timeframes')
    xlabel('Test bin')
    title(['Bins for test ' num2str(C) ' ' num2str(pval)])
end

% decreasing activity
figure('WindowStyle','docked')
for i = 1:4
    subplot(4,3,[1 2]+(i-1)*3)
    [ACinc,AllC,dT] = actionCount(Time{i},Coll{i},4,4,0.05,'Dec');
    stem(dT,ACinc)
    axis([Time{i}(1) Time{i}(end) 0 1])
    xlabel('Time (s)')
    ylabel('Activity-level')
    title([Tag{i} ' rating decreases activity'])
    
    [C,pval,DAct,bins]=sichiSq(ACinc,size(AllC,2),5);
    subplot(4,6,5+(i-1)*6)
    bar([0:size(AllC,2)]/size(AllC,2),DAct)
    legend('Collection','Model','Location','NorthEastOutside')
    ylabel('# timeframes')
    xlabel('Activity-level')
    axis tight
    title('Activity-level distributions')
    
    subplot(4,6,6+(i-1)*6)
    bar(bins)
    legend('Collection','Model','Location','NorthEastOutside')
    ylabel('# timeframes')
    xlabel('Test bin')
    title(['Bins for test ' num2str(C) ' ' num2str(pval)])
end


%% testing joint activities alternating activity


figure('WindowStyle','docked')
for i = 1:4
    subplot(4,5,[1 2]+(i-1)*5)
    [ACinc,AllCinc,dT] = actionCount(Time{i},Coll{i},4,4,0.05,'Inc');
    [ACdec,AllCdec] = actionCount(Time{i},Coll{i},4,4,0.05,'Dec');
    stem(dT,[ACinc -ACdec])
    axis([Time{i}(1) Time{i}(end) -1 1])
    legend Increases Decreases
    xlabel('Time (s)')
    ylabel('Activity-level')
    title([Tag{i} ' rating change activity'])
    
    [C,p,DAct,Bins,v1,v2]=jointChiSq(AllCinc,AllCdec,'Alt');
    bins = [];
    for j = 1:size(Bins,2)
        bins = [bins;squeeze(Bins(:,j,:))];
    end
    
    cmax = max(DAct(DAct>0));
    
    subplot(4,5,3+(i-1)*5)
    imagesc(squeeze(DAct(:,:,1)),[0 cmax])
    axis tight
    xlabel('Total Activity')
    ylabel('Activity Bias')
    title('Joint Expected Activity')
    
    subplot(4,5,4+(i-1)*5)
    imagesc(squeeze(DAct(:,:,2)),[0 cmax])
    axis tight
    xlabel('Total Activity')
    ylabel('Activity Bias')
    title('Joint Actual Activity')
    
    subplot(4,5,5+(i-1)*5)
    bar(bins)
    legend('Collection','Model','Location','NorthEastOutside')
    ylabel('# timeframes')
    xlabel('Test bin')
    title(['Bins for test ' num2str(C) ' ' num2str(p)])
end





%% pairwise inc and dec activity

frameS = 6;

figure('WindowStyle','docked')
for i = 1:2:3
     % inc
    subplot(4,5,[1 2]+(i-1)*5)
    [ACa,AllCa,dT] = actionCount(Time{i},Coll{i},frameS,frameS,0.05,'Inc');
    [ACb,AllCb] = actionCount(Time{i+1},Coll{i+1},frameS,frameS,0.05,'Inc');
    stem(dT,[ACa -ACb])
    axis([Time{i}(1) Time{i}(end) -1 1])
    legend IncreasesA IncreasesB
    xlabel('Time (s)')
    ylabel('Activity-level')
    title([Tag{i} Tag{2} ' rating increases activity'])
    
    [C,p,DAct,Bins,v1,v2]=jointChiSq(AllCa,AllCb,'Ind',3);
    bins = [];
    for j = 1:size(Bins,2)
        bins = [bins;squeeze(Bins(:,j,:))];
    end
    
    cmax = max(DAct(DAct>0));
    
    subplot(4,5,3+(i-1)*5)
    imagesc(squeeze(DAct(:,:,1)),[0 cmax])
    axis tight
    xlabel('Total Activity')
    ylabel('Activity Bias')
    title('Joint Expected Activity')
    
    subplot(4,5,4+(i-1)*5)
    imagesc(squeeze(DAct(:,:,2)),[0 cmax])
    axis tight
    xlabel('Total Activity')
    ylabel('Activity Bias')
    title('Joint Actual Activity')
    
    subplot(4,5,5+(i-1)*5)
    bar(bins)
    legend('Collection','Model','Location','NorthEastOutside')
    ylabel('# timeframes')
    xlabel('Test bin')
    title(['Bins for test ' num2str(C) ' ' num2str(p)])
    
    %Dec
        subplot(4,5,[1 2]+(i)*5)
    [ACa,AllCa,dT] = actionCount(Time{i},Coll{i},frameS,frameS,0.05,'Dec');
    [ACb,AllCb] = actionCount(Time{i+1},Coll{i+1},frameS,frameS,0.05,'Dec');
    stem(dT,[ACa -ACb])
    axis([Time{i}(1) Time{i}(end) -1 1])
    legend IncreasesA IncreasesB
    xlabel('Time (s)')
    ylabel('Activity-level')
    title([Tag{i} Tag{2} ' rating Decreases activity'])
    
    [C,p,DAct,Bins,v1,v2]=jointChiSq(AllCa,AllCb,'Ind',3);
    bins = [];
    for j = 1:size(Bins,2)
        bins = [bins;squeeze(Bins(:,j,:))];
    end
    
    cmax = max(DAct(DAct>0));
    
    subplot(4,5,3+(i)*5)
    imagesc(squeeze(DAct(:,:,1)),[0 cmax])
    axis tight
    xlabel('Total Activity')
    ylabel('Activity Bias')
    title('Joint Expected Activity')
    
    subplot(4,5,4+(i)*5)
    imagesc(squeeze(DAct(:,:,2)),[0 cmax])
    axis tight
    xlabel('Total Activity')
    ylabel('Activity Bias')
    title('Joint Actual Activity')
    
    subplot(4,5,5+(i)*5)
    bar(bins)
    legend('Collection','Model','Location','NorthEastOutside')
    ylabel('# timeframes')
    xlabel('Test bin')
    title(['Bins for test ' num2str(C) ' ' num2str(p)])
end
