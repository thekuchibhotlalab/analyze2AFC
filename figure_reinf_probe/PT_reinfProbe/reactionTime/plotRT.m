%%
%{
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone';

rootPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\';
loadPath = [rootPath filesep 'trialData\'];

allMouse = {};
for i = 1:length(mice)
    load([loadPath mice{i} '.mat']);
    if iscell(task)
        tempObj = wheel2AFC(trialData,'opsPath',[rootPath filesep 'opsFiles' filesep task{i} ],'mouse',mice{i} );
    else
        tempObj = wheel2AFC(trialData,'opsPath',[rootPath filesep 'opsFiles' filesep task ],'mouse',mice{i} );
    end
    tempObj = tempObj.removeMiss(); 
    tempObj = tempObj.getAcc();
    tempObj = tempObj.calculateProbe();
    allMouse{i} = tempObj;
end
%}

clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);

%%
rtReinfBef = {};rtProbe = {}; rtReinf = {}; rtProbeFirstHalf = {}; rtProbeLastHalf = {};
for i = 1:length(mice)
    nDay = unique(allMouse{i}.behav.day);

    for j = 1:length(nDay)
        dayFlag = allMouse{i}.behav.day==nDay(j);
        probeFlag = allMouse{i}.behav.goodProbe==1;
        reinfBefFlag = allMouse{i}.behav.reinfBef==1;

        
        temp = allMouse{i}.behav.reactionTime(~probeFlag & dayFlag,:); 
        tempIdx = round(length(temp)/4):round(length(temp)/4*3); % take the middle 1/2 of the reinf trials of the day
        rtReinf{i}(j) = mean(temp(tempIdx));
        rtReinfBef{i}(j) = mean(allMouse{i}.behav.reactionTime(reinfBefFlag & dayFlag,:));
        
        temp = allMouse{i}.behav.reactionTime(probeFlag & dayFlag,:);
        rtProbe{i}(j) = mean(temp);
        rtProbeFirstHalf{i}(j) = mean(temp(1:round(length(temp)/2))); 
        rtProbeLastHalf{i}(j) = mean(temp(round(length(temp)/2)+1:end)); 
        
        if rtReinf{i}(j)>2.5 || rtReinf{i}(j)<0
            rtReinf{i}(j) = nan; rtProbe{i}(j) = nan; rtReinfBef{i}(j) = nan;
            rtProbeFirstHalf{i}(j) = nan; rtProbeLastHalf{i}(j) = nan;
        end
    end
    
end


figure; [a,b] = fn_sqrtInt(length(mice));
for i = 1:length(mice)
    subplot(a,b,i); hold on; plot(rtReinf{i}); plot(rtReinfBef{i}); hold on; plot(rtProbe{i});
end

rtReinf = fn_cell2matFillNan(rtReinf); rtReinfBef = fn_cell2matFillNan(rtReinfBef); rtProbe = fn_cell2matFillNan(rtProbe);
rtProbeFirstHalf = fn_cell2matFillNan(rtProbeFirstHalf); rtProbeLastHalf = fn_cell2matFillNan(rtProbeLastHalf); 



figure; subplot(3,2,1); p = plotScatter(rtReinf,rtProbe);
title(['reinf RT vs. probe RT; p = ' num2str(p,'%.4f')])
subplot(3,2,2);hold on; p =plotScatter(rtReinfBef,rtProbe);
title(['reinfBef RT vs. probe RT; p = ' num2str(p,'%.4f')])
subplot(3,2,3);hold on; p =plotScatter(rtReinf,rtProbeFirstHalf);
title(['reinf RT vs. probeFirstHalf RT; p = ' num2str(p,'%.4f')]);
subplot(3,2,4);hold on; p =plotScatter(rtReinfBef,rtProbeFirstHalf);
title(['reinfBef RT vs. probeFirstHalf RT; p = ' num2str(p,'%.4f')]);
subplot(3,2,5);hold on; p = plotScatter(rtReinf,rtProbeLastHalf);
title(['reinf RT vs. probeLastHalf RT; p = ' num2str(p,'%.4f')]);
subplot(3,2,6);hold on; p = plotScatter(rtReinfBef,rtProbeLastHalf);
title(['reinfBef RT vs. probeLastHalf RT; p = ' num2str(p,'%.4f')]);

figure; subplot(3,2,1);hold on;  plot(nanmean(rtReinf,1)); plot(nanmean(rtProbe,1)); title('reinf vs. probe')
subplot(3,2,2);hold on; plot(nanmean(rtReinfBef,1)); plot(nanmean(rtProbe,1)); title('reinfBef vs. probe')
subplot(3,2,3);hold on; plot(nanmean(rtReinf,1)); plot(nanmean(rtProbeFirstHalf,1)); title('reinf vs. probeFirstHalf')
subplot(3,2,4);hold on; plot(nanmean(rtReinfBef,1)); plot(nanmean(rtProbeFirstHalf,1)); title('reinfBef vs. probeFirstHalf')
subplot(3,2,5);hold on; plot(nanmean(rtReinf,1)); plot(nanmean(rtProbeLastHalf,1)); title('reinf vs. probeLastHalf')
subplot(3,2,6);hold on; plot(nanmean(rtReinfBef,1)); plot(nanmean(rtProbeLastHalf,1)); title('reinfBef vs. probeLastHalf')


%%
for i = 1:length(mice)
    figure; 
    nDay = unique(allMouse{i}.behav.day);
    [a,b] = fn_sqrtInt(length(nDay));
    for j = 1:length(nDay)
        dayFlag = allMouse{i}.behav.day==nDay(j);
        reReinf{j} = smoothdata(allMouse{i}.behav.reactionTime(dayFlag,:),'movmean',30);
        subplot_tight(a,b,j); plot(reReinf{j});     
    end
end

%% look at correlation of RT and bias, combining both left and right
%mouse = 1; day = 1; nTrials = 27:117;
%mouse = 11; day = 4; nTrials = 65:264;
mouse = 11; day = 3; nTrials = 43:242;
%mouse = 2; day = 8; nTrials = 1:291;
%mouse = 6; day = 4; nTrials = 36:235;
%mouse = 13; day = 8; nTrials = 18:257;
dayFlag = (allMouse{mouse}.behav.day==day);
tempRT = allMouse{mouse}.behav.reactionTime(dayFlag,:); tempRT = tempRT(nTrials);
tempBias = allMouse{mouse}.behav.bias(dayFlag,:); tempBias = tempBias(nTrials);
%tempBias = abs(allMouse{mouse}.behav.bias(dayFlag,:)); tempBias = tempBias(nTrials);


nBins = 3;
ymax = (max(tempBias)/0.1)*0.1; ymin = (min(tempBias)/0.1)*0.1; 
ybins = linspace(ymin,ymax,nBins+1);
figure; subplot(4,1,1); plot(tempRT);
subplot(4,1,2); plot(tempBias)
subplot(2,2,3)
scatter(tempRT,tempBias); ylim([ymin ymax]);
subplot(2,2,4); hold on;
[~,~,binIdx] = histcounts(tempBias,ybins); tempRTbin = cell(1,nBins);
for i = 1:nBins; tempRTbin{i} = tempRT(binIdx==i); cdfplot(tempRTbin{i}); end 

%% look at correlation of RT and bias, separating both left and right
%mouse = 1; day = 1; nTrials = 27:117;
mouse = 11; day = 4; nTrials = 65:264;
%mouse = 11; day = 3; nTrials = 43:242;
%mouse = 2; day = 8; nTrials = 1:291;
%mouse = 6; day = 4; nTrials = 36:235;
%mouse = 13; day = 8; nTrials = 18:257;
dayFlag = (allMouse{mouse}.behav.day==day);
tempRT = allMouse{mouse}.behav.reactionTime(dayFlag,:); tempRT = tempRT(nTrials);
tempBias = allMouse{mouse}.behav.bias(dayFlag,:); tempBias = tempBias(nTrials);
tempAction = allMouse{mouse}.behav.action(dayFlag,:); tempAction = tempAction(nTrials);
%tempBias = abs(allMouse{mouse}.behav.bias(dayFlag,:)); tempBias = tempBias(nTrials);
plotRT_LR(tempRT,tempBias,tempAction);

%% look at all days (all trials) for each animals
mouse = 13;
nDay = unique(allMouse{mouse}.behav.day);

for i = 1:length(nDay)
    dayFlag = (allMouse{mouse}.behav.day==nDay(i));
    tempRT = allMouse{mouse}.behav.reactionTime(dayFlag,:); %tempRT = tempRT - smoothdata(tempRT,'movmean',200);
    tempBias = allMouse{mouse}.behav.bias(dayFlag,:); 
    tempAction = allMouse{mouse}.behav.action(dayFlag,:);
    plotRT_LR(tempRT,tempBias,tempAction,[0 2.5]);
end

%% look at the wheel traj for all days of each animal
mouse = 9;
nDay = unique(allMouse{mouse}.behav.day);

for i = 1:length(nDay)
    dayFlag = (allMouse{mouse}.behav.day==nDay(i));
    tempWheel = allMouse{mouse}.behavVar.wheelPos_aligned(dayFlag,:);
    tempRT = allMouse{mouse}.behav.reactionTime(dayFlag,:);
    tempBias = allMouse{mouse}.behav.bias(dayFlag,:); 
    tempAction = allMouse{mouse}.behav.action(dayFlag,:);
    plotWheel_LR(tempRT, tempWheel,tempBias,tempAction)
end
%% look at the wheel traj only for RT>1sec
mouse = 9;
highRT = allMouse{mouse}.behav.reactionTime > 1;
tempWheel = allMouse{mouse}.behavVar.wheelPos_aligned(highRT,:);
figure; plot(tempWheel(1:10,:)')
%% characterize wheel trajectory 
mouse = 9;
[wheelPos, badFlag] = removeBadTrials(allMouse{mouse}.behavVar.wheelPos_aligned);
nTrials = size(wheelPos,1);
behav = allMouse{mouse}.behav(~badFlag,:);
tempTimeDownSampleIdx = 2001:50:4501;
tempWheelDownSample = wheelPos(:,tempTimeDownSampleIdx);
tempWheelChange = diff(tempWheelDownSample,1,2);
wheelOneBlockIdx = false(nTrials,1); wheelMultiBlockIdx = false(nTrials,1); tempOnsetIdx = zeros(nTrials,1);
for i = 1:nTrials
    % find one-block or multiple block movements
    [onIdx, offIdx,blockIdx] = fn_getBlockOnOff(tempWheelChange(i,:)~=0);
    if length(onIdx)>1; wheelMultiBlockIdx(i) = 1;
    else; wheelOneBlockIdx(i) = 1; end      
    % align and do clustering
    tempOnsetIdx(i) = find(abs(tempWheelDownSample(i,:))~=0,1); 
end
[tempWheelDownSample_aligned,maxAlignPoint] = fn_align2idx(tempWheelDownSample, tempOnsetIdx,'fill','startEndValue');
figure;
subplot(4,1,1); hold on;  plot([2000 2000],[-40 40],'Color',[0.8 0.8 0.8],'LineWidth',2); 
plot(wheelPos([2 4],:)','Color',[0.2 0.2 0.2],'LineWidth',1); xlim([1500 4500]); ylim([-40 40]); 
xticks([2000 3000 4000]); xticklabels({'0','1','2'})
yticks([-35 0 35]); yticklabels({'R','0','L'})
subplot(4,1,2); hold on;  plot([2000 2000],[-40 40],'Color',[0.8 0.8 0.8],'LineWidth',2); 
plot(wheelPos(1:30,:)','Color',[0.2 0.2 0.2],'LineWidth',1); xlim([1500 4500]); ylim([-40 40])
xticks([2000 3000 4000]); xticklabels({'0','1','2'})
yticks([-35 0 35]); yticklabels({'R','0','L'})
subplot(4,1,3); hold on;  plot([2000 2000],[-40 40],'Color',[0.8 0.8 0.8],'LineWidth',2); 
plot(wheelPos(wheelOneBlockIdx(1:30),:)','Color',matlabColors(1),'LineWidth',1); xlim([1500 4500]); ylim([-40 40]); xticks([])
yticks([-35 0 35]); yticklabels({'R','0','L'})
subplot(4,1,4); hold on; plot([2000 2000],[-40 40],'Color',[0.8 0.8 0.8],'LineWidth',2); 
plot(wheelPos(wheelMultiBlockIdx(1:30),:)','Color',matlabColors(2),'LineWidth',1); xlim([1500 4500]);  ylim([-40 40])
xticks([2000 3000 4000]); xticklabels({'0','1','2'}); yticks([-35 0 35]); yticklabels({'R','0','L'})

figure;  subplot(2,1,1);plot(tempWheelDownSample_aligned(wheelOneBlockIdx(1:100),:)','Color',matlabColors(1),'LineWidth',1)
xlim([40  100]);xticks([]);ylim([-40 40]);yticks([-35 0 35]); yticklabels({'R','0','L'})
subplot(2,1,2);  plot(tempWheelDownSample_aligned(wheelMultiBlockIdx(1:30),:)','Color',matlabColors(2),'LineWidth',1)
tempWheelDownSample_aligned = tempWheelDownSample_aligned(:,maxAlignPoint-1:end);
xlim([40 100]);xticks([]);ylim([-40 40]);yticks([-35 0 35]); yticklabels({'R','0','L'})

% change the sign of aligned wheel traj to match 
tempChangeFlag = tempWheelDownSample_aligned(:,end)<0; 
tempWheelDownSample_aligned(tempChangeFlag,:) = -tempWheelDownSample_aligned(tempChangeFlag,:);

plotWheelRT_cluster(behav,{wheelOneBlockIdx, wheelMultiBlockIdx},{'Fast', 'Slow'})
%% do a t-sne
Y = tsne(tempWheelDownSample_aligned);
%% do k means

k = 5; idx = kmeans(abs(tempWheelDownSample_aligned),k);
%% visualization of kmeans 
figure;
colors = parula; tempIdx = round(linspace(1, size(colors,1),k+1));
colors = colors(tempIdx(1:end-1),:);
for i = 1:k;  subplot(k,1,i);
    tempIdx = find(idx==i); tempPlotIdx = randperm(length(tempIdx));
    plot(tempWheelDownSample_aligned(tempIdx(tempPlotIdx(1:110)),:)','Color',colors(i,:),'LineWidth',0.2);
xticks([]);ylim([-40 40]);yticks([-35 0 35]); yticklabels({'R','0','L'})
end

sampleColor = nan(size(tempWheelDownSample_aligned,1),3);
for i = 1:k; sampleColor(idx==i,:) = repmat(colors(i,:),sum(idx==i),1); end 
figure; scatter(Y(:,1),Y(:,2),30,sampleColor,'filled'); xticks([]); yticks([])
% look at cluster distribution across learning trials
nTrials = size(wheelPos,1);
tempBinSize = 600; tempBin = 0: tempBinSize: nTrials; tempBin(end) = nTrials;
tempCounts = zeros(k,length(tempBin)-1); for i = 1:k; [tempCounts(i,:)] = histcounts(find(idx==i),tempBin); end 
tempCounts = fn_normalizeBySum(tempCounts,2);
figure; plot(tempCounts')


%% Group the clusters
clusterGroup = {4,1,[2 3 5]};clusterGroupLabel = {'Fast','Mid','Slow'}; clusterGroupFlag = {}; 
for i = 1: length(clusterGroup); clusterGroupFlag{i} = sum(idx==clusterGroup{i},2)>0; end
% Bias and RT analysis for different type of wheel trajectory
plotWheelRT_cluster(behav,clusterGroupFlag,clusterGroupLabel)

%% align trajectory to movement onset for one-block movements
tempWheelOneBlock = wheelPos(wheelOneBlockIdx,:);tempOnsetIdx = [];
for i = 1:size(tempWheelOneBlock,1)
    tempOnsetIdx(i) = find(abs(tempWheelOneBlock(i,2001:end))~=0,1)+2000; 
end
[tempMat,~] = fn_align2idx(tempWheelOneBlock, tempOnsetIdx);
figure; plot(tempMat')

%% characterize wheel trajectory 
mice = 8:13;nTrialMax = 5000;
wheelOneBlockIdx = nan(nTrialMax,length(mice)); wheelMultiBlockIdx = nan(nTrialMax,length(mice)); tempOnsetIdx = nan(nTrialMax,length(mice));
for j = 1:length(mice)
    mouse = mice(j);
    [wheelPos, badFlag] = removeBadTrials(allMouse{mouse}.behavVar.wheelPos_aligned);
    nTrials = size(wheelPos,1);
    behav = allMouse{mouse}.behav(~badFlag,:);
    tempTimeDownSampleIdx = 2001:50:4501;
    tempWheelDownSample = wheelPos(:,tempTimeDownSampleIdx);
    tempWheelChange = diff(tempWheelDownSample,1,2);
    
    for i = 1:nTrials
        % find one-block or multiple block movements
        [onIdx, offIdx,blockIdx] = fn_getBlockOnOff(tempWheelChange(i,:)~=0);
        if length(onIdx)>1; wheelMultiBlockIdx(i,j) = 1; wheelOneBlockIdx(i,j) = 0;
        else; wheelOneBlockIdx(i,j) = 1; wheelMultiBlockIdx(i,j) = 0; end      
        % align and do clustering
        tempOnsetIdx(i,j) = find(abs(tempWheelDownSample(i,:))~=0,1); 
    end
end
figure; hold on; 
for i = 1:length(mice); plot(smoothdata(wheelOneBlockIdx,'movmean',200, 'includenan' ),'Color',[0.8 0.8 0.8],'LineWidth',1); end
plot(nanmean(smoothdata(wheelOneBlockIdx,'movmean',200, 'includenan' ),2),'Color',[0.2 0.2 0.2],'LineWidth',2)
xlim([0 2500])

%% all functions
function [tempWheel, badFlag] = removeBadTrials(tempWheel)
    totalTrial = size(tempWheel,1);
    badFlag = tempWheel(:,end)==0;
    tempWheel = tempWheel(~badFlag,:);
    disp(['Bad Trials = ' int2str(sum(badFlag)) ' out of ' int2str(totalTrial)])
end



function plotRT_LR(tempRT,tempBias,tempAction,xlimm)
if nargin == 3; xlimm = [0 2.5]; end 
tempRTL = tempRT(tempAction==1); tempRTR = tempRT(tempAction==2);
tempBiasL = tempBias(tempAction==1); tempBiasR = tempBias(tempAction==2);

nBins = 3;
ymax = max(tempBiasL); ymin = min(tempBiasL); 
ybins = linspace(ymin,ymax,nBins+1); 
figure; subplot(6,1,1); plot(tempRT);
subplot(6,1,2); plot(tempBias)
subplot(3,2,3)
scatter(tempRTL,tempBiasL); ylim([ymin ymax]); xlim(xlimm)
subplot(3,2,4); hold on;
[~,~,binIdx] = histcounts(tempBiasL,ybins); tempRTbin = cell(1,nBins);
for i = 1:nBins; tempRTbin{i} = tempRTL(binIdx==i); cdfplot(tempRTbin{i}); end; xlim(xlimm)

nBins = 3;
ymax = (max(tempBiasR)/0.1)*0.1; ymin = (min(tempBiasR)/0.1)*0.1; 
ybins = linspace(ymin,ymax,nBins+1);
subplot(3,2,5)
scatter(tempRTR,tempBiasR); ylim([ymin ymax]); xlim([0 2.5])
subplot(3,2,6); hold on;
[~,~,binIdx] = histcounts(tempBiasR,ybins); tempRTbin = cell(1,nBins);
for i = 1:nBins; tempRTbin{i} = tempRTR(binIdx==i); cdfplot(tempRTbin{i}); end; xlim(xlimm)    
end

% plot wheel traja nd RT given a set of clusters
function plotWheelRT_cluster(behav,clusterGroupFlag,clusterGroupLabel)
biasThreshold = 0.2;
highBiasFlag = behav.bias>=biasThreshold | behav.bias<=-biasThreshold;
lowBiasFlag = ~highBiasFlag;
    
figure; subplot(4,2,1); hold on; for i = 1:length(clusterGroupFlag); plot(smoothdata(clusterGroupFlag{i},'movmean',300),'LineWidth',2);end
legend(clusterGroupLabel); ylim([0 1])


subplot(4,2,3); hold on;
plot(smoothdata(clusterGroupFlag{1}&highBiasFlag,'movmean',300)./smoothdata(clusterGroupFlag{1},'movmean',300),'LineWidth',2)
plot(smoothdata(clusterGroupFlag{2}&highBiasFlag,'movmean',300)./smoothdata(clusterGroupFlag{2},'movmean',300),'LineWidth',2)


subplot(4,2,2); hold on; 
for i = 1:length(clusterGroupFlag)
    tempRT = nan(size(behav.reactionTime)); tempRT(clusterGroupFlag{i}) = behav.reactionTime(clusterGroupFlag{i});
    plot(smoothdata(tempRT,'movmean',150),'LineWidth',2);
end
legend({'High Bias','Low Bias'});ylim([0 2]); yticks(0:0.5:2)


subplot(4,2,4); hold on; plot(smoothdata(behav.actionRate,'movmean',1));


tempBiasCluster = cellfun(@(x)behav.bias(x),clusterGroupFlag,'UniformOutput',false);
subplot(2,2,3);hold on; cellfun(@cdfplot,tempBiasCluster);legend(clusterGroupLabel)
tempRTCluster = cellfun(@(x)behav.reactionTime(x),clusterGroupFlag,'UniformOutput',false);
subplot(2,2,4); hold on; h = cellfun(@cdfplot,tempRTCluster,'UniformOutput',false);legend(clusterGroupLabel)
cellfun(@(x)(set(x,'LineWidth',2)),h); title('RT Distribution'); xlabel('Reaction Time'); ylabel('Cumulative Prop')

choiceL = behav.action==1; choiceR = behav.action==2;
for i = 1:length(clusterGroupFlag)
    figure;
    highBiasL = behav.bias>=biasThreshold & choiceL & clusterGroupFlag{i}; antiBiasL = behav.bias<=-biasThreshold & choiceL & clusterGroupFlag{i};
    noBiasL = (behav.bias>-biasThreshold & behav.bias<biasThreshold)& choiceL & clusterGroupFlag{i};
    highBiasR = behav.bias<=-biasThreshold & choiceR & clusterGroupFlag{i}; antiBiasR = behav.bias>=biasThreshold & choiceR & clusterGroupFlag{i};
    noBiasR = (behav.bias>-biasThreshold & behav.bias<biasThreshold)& choiceR & clusterGroupFlag{i};
    subplot(2,1,1); hold on; 
    cdfplot(behav.reactionTime(highBiasL)); cdfplot(behav.reactionTime(antiBiasL)); cdfplot(behav.reactionTime(noBiasL))
    [~,p] = ttest2(behav.reactionTime(highBiasL),behav.reactionTime(antiBiasL));
    title(['Wheel Group ' int2str(i) '. Left choice, p = ' num2str(p,'%.4f')])
    legend({['L act L bias, mean = ' num2str(mean(behav.reactionTime(highBiasL)),'%.2f')],...
        ['L act R bias, mean = ' num2str(mean(behav.reactionTime(antiBiasL)),'%.2f')],...
        ['L act no bias, mean = ' num2str(mean(behav.reactionTime(noBiasL)),'%.2f')]},'Location','Best')

    subplot(2,1,2); hold on; 
    cdfplot(behav.reactionTime(highBiasR)); cdfplot(behav.reactionTime(antiBiasR)); cdfplot(behav.reactionTime(noBiasR))
    [~,p] = ttest2(behav.reactionTime(highBiasR),behav.reactionTime(antiBiasR));
    title(['Wheel Group ' int2str(i) '. Right choice, p = ' num2str(p,'%.4f')])
    legend({['R act R bias, mean = ' num2str(mean(behav.reactionTime(highBiasR)),'%.2f')],...
        ['R act L bias, mean = ' num2str(mean(behav.reactionTime(antiBiasR)),'%.2f')],...
        ['R act no bias, mean = ' num2str(mean(behav.reactionTime(noBiasR)),'%.2f')]},'Location','Best')
end
    
figure;
for i = 1:length(clusterGroupFlag)
    subplot(length(clusterGroupFlag),1,i); hold on;
    tempHighFlag = highBiasFlag & clusterGroupFlag{i} & choiceL; tempLowFlag = lowBiasFlag & clusterGroupFlag{i} & choiceL;
    highBiasRT = behav.reactionTime; highBiasRT(~tempHighFlag) = nan;
    lowBiasRT = behav.reactionTime; lowBiasRT(~tempLowFlag) = nan;
    plot(smoothdata(highBiasRT,'movmean',500));plot(smoothdata(lowBiasRT,'movmean',500));
end

end

% plot wheel trajectory
function plotWheel_LR(tempRT, tempWheel,tempBias,tempAction)
wheelAxis = 1001:4501;
tempWheel = tempWheel(:,wheelAxis);

figure; subplot(6,1,1); plot(tempRT); xlim([1 length(tempRT)])
subplot(6,1,2); plot(tempBias); xlim([1 length(tempRT)])

nBins = 3;
ymax = max(tempBias); ymin = min(tempBias); 
ybins = linspace(ymin,ymax,nBins+1); 

[~,~,binIdx] = histcounts(tempBias,ybins); tempWheelBin = cell(1,nBins);
for i = 1:nBins
    
    tempWheelBin{i} = tempWheel(binIdx==i,:); 
    tempYMax = max(tempWheelBin{i}(:,2001:end)); tempYMax = max(tempYMax(:));
    tempWheelBin{i} = tempWheelBin{i}/tempYMax;
    tempL = tempAction(binIdx==i); tempL = tempL==1;
    tempR = tempAction(binIdx==i); tempR = tempR==2;
    
    subplot(6,1,3)
    fn_plotMeanErrorbar(1:length(wheelAxis),tempWheelBin{i}(tempL,:),matlabColors(i),{'Color',matlabColors(i)},{'faceAlpha',0.3,'LineStyle','none'});
    fn_plotMeanErrorbar(1:length(wheelAxis),tempWheelBin{i}(tempR,:),matlabColors(i),{'Color',matlabColors(i)},{'faceAlpha',0.3,'LineStyle','none'});
    xlim([901 1501]); xticks([1 1001 1501 2001 2501 3001]); xticklabels({'-1','0','0.5','1','1.5','2'}); ylim([-1.2 1.2]);
    title(['Bias Bin From ' num2str(ybins(i),'%.2f') ' to ' num2str(ybins(i+1),'%.2f')])
    
    subplot(6,1,7-i)
    plot(tempWheelBin{i}')
    xlim([901 1501]); xticks([1 1001 1501 2001 2501 3001]); xticklabels({'-1','0','0.5','1','1.5','2'}); ylim([-1.2 1.2]);
    title(['Bias Bin From ' num2str(ybins(i),'%.2f') ' to ' num2str(ybins(i+1),'%.2f')])
end
  
end

function p = plotScatter(x,y)
hold on; plot([0 1.5],[0 1.5],'Color',[0.8 0.8 0.8])
scatter(x(:),y(:),50,matlabColors(1),'filled'); xlim([0 1.5]); ylim([0 1.5])
[~,p] = ttest(x(:),y(:)); 
    
end