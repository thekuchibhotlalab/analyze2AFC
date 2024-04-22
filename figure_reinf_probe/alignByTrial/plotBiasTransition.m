function plotBiasTransition(mice,selectProtocol)

if nargin == 1; selectProtocol = {'puretone'}; end

rootPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\';
loadPath = [rootPath filesep 'trialData\'];
figPath = [rootPath filesep 'figure\' filesep 'learningCurveTrialAll' filesep]; mkdir(figPath);
trialData = cellfun(@(x)fn_selectProtocol(loadPath, x, selectProtocol),mice,'UniformOutput',false);

trialData = cellfun(@fn_selectAnimalDay,trialData,mice,'UniformOutput',false);

%[~, accuracy,~,~,~,~,~,~,~,~,~,~,probeData,~, ~,~,probeTrialNum,~,~]
outmat = cellfun(@(x)alignDays(x),trialData,'UniformOutput',false);


[L2R,R2L,midL2R,midR2L] = cellfun(@computeTrials,outmat,mice,'UniformOutput',false);

% PLOT REACTION TIME IN BLOCK
meanRT = cellfun(@plotTransActionResponseTime,outmat,L2R,R2L, mice,'UniformOutput',false);
meanRT = fn_cell2mat(meanRT,1); meanRTplot = [mean(meanRT(:,1:2),2)  mean(meanRT(:,3:4),2)];
figure; hold on;
bar(mean(meanRTplot,1),'EdgeColor','None','FaceColor',matlabColors(2,0.9))
for i = 1:size(meanRTplot,1)
    f = plot([1 2],[meanRTplot(i,1) meanRTplot(i,2)],'Color',[0.6 0.6 0.6],'Marker','.','MarkerSize',15,...
        'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',0.5);
end
xticks([1 2]); xlim([0 3]); xticklabels({'Biased Block', 'Unbiased Block'});xtickangle(25);
ylabel('Reaction Time (ms)')

% PLOT BLOCK TRANSITION FREQ
L2R = fn_cell2mat(L2R,1); R2L = fn_cell2mat(R2L,1);
midL2R = fn_cell2mat(midL2R,1); midR2L = fn_cell2mat(midR2L,1);
figure; histogram([midL2R;midR2L],50)
xlim([0 2500]);xlabel('Trials'); ylabel('Frequency')
title('Frequency of Bias Transitions')

figure;
scatter((L2R(:,1) + L2R(:,2))/2, L2R(:,2) - L2R(:,1))
xlabel('Trials'); ylabel('Bias block len'); title('L2R')

figure;
scatter((R2L(:,1) + R2L(:,2))/2, R2L(:,2) - R2L(:,1))
xlabel('Trials'); ylabel('Bias block len'); title('R2L')





blockIdxL = fn_cell2matFillNan(cellfun(@(x)(fn_idx2logical(fn_catIdx([x.biasBlockL_start x.biasBlockL_end]),length(x.stimulus))),...
    outmat,'UniformOutput',false));

blockIdxR = fn_cell2matFillNan(cellfun(@(x)(fn_idx2logical(fn_catIdx([x.biasBlockR_start x.biasBlockR_end]),length(x.stimulus))),...
    outmat,'UniformOutput',false));
blockIdx = blockIdxL + blockIdxR; blockIdx(blockIdx>1) = 1;
bins = 0:400:2400;
tempSum = [];
for i = 1:length(bins)-1
    tempSum(:,i) = sum(blockIdx((bins(i)+1):bins(i+1),:),1);
end
%tempSum(6,:) = [];
figure; hold on

fn_plotFillErrorbar(bins(1:end-1)+bins(2)/2,nanmean(tempSum,1), nanstd(tempSum,0,1)./ sqrt(size(tempSum,1)),...
    matlabColors(1),'FaceAlpha',0.2,'LineStyle','None')
plot(bins(1:end-1)+bins(2)/2,nanmean(tempSum,1),'Color',matlabColors(1))

midL = cellfun(@(x)(round((x.biasBlockL_start+x.biasBlockL_end)/2)),outmat,'UniformOutput',false);
midR = cellfun(@(x)(round((x.biasBlockR_start+x.biasBlockR_end)/2)),outmat,'UniformOutput',false);
lenL = cellfun(@(x)(x.biasBlockL_len),outmat,'UniformOutput',false);
lenR = cellfun(@(x)(x.biasBlockR_len),outmat,'UniformOutput',false);

%bins = 0:400:2800;
blockLen = nan(length(midL),length(bins)-1);
for i = 1:length(midL)
    [~,~,whichBinL] = histcounts(midL{i},bins); [~,~,whichBinR] = histcounts(midR{i},bins);
    for j = 1:length(bins)-1
        tempL = lenL{i}(whichBinL == j); tempR = lenR{i}(whichBinR == j);
        blockLen(i,j) = mean([tempL;tempR]);
    end
end
blockLen(6,:) = [];
figure; hold on;
fn_plotFillErrorbar(bins(1:end-1)+bins(2)/2,nanmean(blockLen,1), nanstd(blockLen,0,1)./ sqrt(size(blockLen,1)),...
    matlabColors(1),'FaceAlpha',0.2,'LineStyle','None')
plot(bins(1:end-1)+bins(2)/2,nanmean(blockLen,1),'Color',matlabColors(1))



figure; subplot(2,1,1); hold on;

tempSum = tempSum ./ mean(mean(tempSum(:,1:2))); blockLen = blockLen ./ mean(mean(blockLen(:,1:2)));
scatter(bins(1:end-1)+bins(2)/2,nanmean(tempSum,1),20,matlabColors(1),'filled')


scatter(bins(1:end-1)+bins(2)/2,nanmean(blockLen,1),20,matlabColors(6),'filled')
f = lsline;
set(f(1), 'Color', matlabColors(6),'LineWidth',2); set(f(2), 'Color', matlabColors(1),'LineWidth',2)
legend(f,{'Bias Block Len','Bias Block Total'})
ylabel('Normalized Unit'); xlabel('Trials');ylim([0.2 1.4]); yticks([0.2:0.4:1.4]); xlim([0 2500])

subplot(2,1,2); hold on;
midTrans = [midL2R;midR2L]; midTransFreq = histcounts(midTrans,bins)./sum(~isnan(tempSum),1);
scatter(bins(1:end-1)+bins(2)/2,midTransFreq./mean(midTransFreq(1:2)),20,matlabColors(2),'filled');
f = lsline;

set(f(1), 'Color', matlabColors(2),'LineWidth',2)
legend(f,{'Transitions'})
ylabel('Normalized Unit'); xlabel('Trials');ylim([0.2 1.4]); yticks([0.2:0.4:1.4]); xlim([0 2500])

end

function [blockTrialL2R,blockTrialR2L,midPointL2R,midPointR2L] = computeTrials(outmat,mouse)

action = outmat.action; action(action==2) = 0; actionBin = 30;
action = smoothdata(action,'movmean',actionBin);

stimulus = outmat.stimulus; stimulus(stimulus==2) = 0; 
stimulus = smoothdata(stimulus,'movmean',actionBin);

tempL_start = outmat.biasBlockL_start; tempR_end = outmat.biasBlockR_end;
tempL_end = outmat.biasBlockL_end; tempR_start = outmat.biasBlockR_start;

% Plot L2R transitions
tempDiff = repmat(tempR_start',length(tempL_end),1)-repmat(tempL_end,1,length(tempR_start)); 
for i = 1:length(tempL_end)
    tempIdx = find(tempDiff(i,:) > 0); 
    if ~isempty(tempIdx)
        tempIdx = tempIdx(1);
        blockDiffL2R(i) = tempDiff(i,tempIdx);
        blockTrialL2R(i,1) = tempL_start(i); blockTrialL2R(i,2) = tempL_end(i); 
        blockTrialL2R(i,3) = tempR_start(tempIdx); blockTrialL2R(i,4) = tempR_end(tempIdx);
    end
end
transTrialThreshold = 30; thresholdFlag = (blockTrialL2R(:,3)-blockTrialL2R(:,2)) <= transTrialThreshold;
blockTrialL2R = blockTrialL2R(thresholdFlag,:);

plotTransActionStim(blockTrialL2R,action,stimulus,[mouse ' L2R ']);


% Plot R2L transitions
tempDiff = repmat(tempL_start',length(tempR_end),1)-repmat(tempR_end,1,length(tempL_start)); 
for i = 1:length(tempR_end)
    tempIdx = find(tempDiff(i,:) > 0); 
    if ~isempty(tempIdx)
        tempIdx = tempIdx(1);
        blockDiffR2L(i) = tempDiff(i,tempIdx);
        blockTrialR2L(i,1) = tempR_start(i); blockTrialR2L(i,2) = tempR_end(i); 
        blockTrialR2L(i,3) = tempL_start(tempIdx); blockTrialR2L(i,4) = tempL_end(tempIdx);
    end
end
transTrialThreshold = 30; thresholdFlag = (blockTrialR2L(:,3)-blockTrialR2L(:,2)) <= transTrialThreshold;
blockTrialR2L = blockTrialR2L(thresholdFlag,:);

plotTransActionStim(blockTrialR2L,action,stimulus,[mouse ' R2L ']);

midPointL2R = round((blockTrialL2R(:,2)+blockTrialL2R(:,3))/2);
midPointR2L = round((blockTrialR2L(:,2)+blockTrialR2L(:,3))/2);

end

function meanRT = plotTransActionResponseTime(trialData,L2R,R2L, mouse)

responseTime = trialData.responseTime;
responseType = trialData.responseType;

temp = 0;
startL = trialData.biasBlockL_start-temp; startL(startL<1)=1;
endL = trialData.biasBlockL_end+temp; endL(endL>length(responseTime))=length(responseTime);
idxL = fn_catIdx([startL endL]); idxL = fn_idx2logical(idxL,length(responseTime));

startR = trialData.biasBlockR_start-temp; startR(startR<1)=1;
endR = trialData.biasBlockR_end+temp; endR(endR>length(responseTime))=length(responseTime);
idxR = fn_catIdx([startR endR]); idxR = fn_idx2logical(idxR,length(responseTime));

if ~isempty(L2R)    
    [inBlockL2R, inTransL2R] = getBlockIdx(L2R,length(responseTime)); 
    [inBlockR2L, inTransR2L] = getBlockIdx(R2L,length(responseTime)); 
    
    inBlock = inBlockL2R | inBlockR2L;
    inTrans = inTransL2R | inTransR2L;
    
    inBlockRT_corr = responseTime(inBlock & responseType==1);
    inBlockRT_incorr = responseTime(inBlock & responseType==2);

    inTransRT_corr = responseTime(inTrans & responseType==1);
    inTransRT_incorr = responseTime(inTrans & responseType==2);

    outBlock = ~(inBlock & inTrans);
    outBlockRT_corr = responseTime(outBlock & responseType==1);
    outBlockRT_incorr = responseTime(outBlock & responseType==2);

    figure; 
    meanRT = [mean(inBlockRT_corr),mean(inBlockRT_incorr),mean(outBlockRT_corr),...
        mean(outBlockRT_incorr),mean(inTransRT_corr),mean(inTransRT_incorr)];
    bar(meanRT);
    temp = [length(inBlockRT_corr),length(inBlockRT_incorr),length(outBlockRT_corr),...
        length(outBlockRT_incorr),length(inTransRT_corr),length(inTransRT_incorr)];
    xticklabels(strsplit(int2str(temp))); title(mouse)
    
else
    inTrans = false(length(responseTime),1);
    meanRT = [];
end
%{
inBlockRT_corr = responseTime((idxL | idxR) & ~(inTrans) & responseType==1);
inBlockRT_incorr = responseTime((idxL | idxR) & ~(inTrans) & responseType==2);

outBlockRT_corr = responseTime(~(idxL | idxR) & ~(inTrans) & responseType==1);
outBlockRT_incorr = responseTime(~(idxL | idxR) & ~(inTrans) & responseType==2);

figure; 
bar([mean(inBlockRT_corr), mean(inBlockRT_incorr),mean(outBlockRT_corr),mean(outBlockRT_incorr)])
%}

%disp('done')

end

function [inBlock, inTrans] = getBlockIdx(L2R,len)
    inBlock = []; inTrans = [];
    for i = 1:size(L2R)
        inBlock = [inBlock L2R(i,1):L2R(i,2) L2R(i,3):L2R(i,4)];
        inTrans = [inTrans L2R(i,2):L2R(i,3)];
    end
     inBlock = fn_idx2logical(inBlock,len); 
     inTrans = fn_idx2logical(inTrans,len);
end

function plotTransActionStim(blockTrialL2R,action,stimulus,mouse)

%{
figure; subplot(1,2,1); hold on
for i = 1:size(blockTrialL2R,1)
    tempTrialIdx = blockTrialL2R(i,1):blockTrialL2R(i,4);
    midPoint = (blockTrialL2R(i,2)+blockTrialL2R(i,3))/2; 
    plot(tempTrialIdx - midPoint, -i*0.5 + action(tempTrialIdx))
    
end
xlim([-30 30])

subplot(1,2,2); hold on
for i = 1:size(blockTrialL2R,1)
    tempTrialIdx = blockTrialL2R(i,1):blockTrialL2R(i,4);
    midPoint = (blockTrialL2R(i,2)+blockTrialL2R(i,3))/2;
    
    plot(tempTrialIdx - midPoint, -i*0.5 + stimulus(tempTrialIdx))
    
end
xlim([-30 30])
%}

tempTrialIdxMax = [];tempTrialIdxMin = [];
for i = 1:size(blockTrialL2R,1)
    midPoint(i) = round((blockTrialL2R(i,2)+blockTrialL2R(i,3))/2);
    tempTrialIdxMax(i) = blockTrialL2R(i,4) - midPoint(i);
    tempTrialIdxMin(i) = midPoint(i) - blockTrialL2R(i,1);
end
if ~isempty(tempTrialIdxMax) && ~isempty(tempTrialIdxMin)
    tempTrialIdxMax = max(tempTrialIdxMax); tempTrialIdxMin = max(tempTrialIdxMin);
else
    tempTrialIdxMax = nan; tempTrialIdxMin = nan;
end

if ~isnan(tempTrialIdxMax)

    actionTransPlot = ones(size(blockTrialL2R,1),tempTrialIdxMax+tempTrialIdxMin+1)*0.5;
    stimulusTransPlot = ones(size(blockTrialL2R,1),tempTrialIdxMax+tempTrialIdxMin+1)*0.5;
    
    actionTransPlotNan = nan(size(blockTrialL2R,1),tempTrialIdxMax+tempTrialIdxMin+1)*0.5;
    stimulusTransPlotNan = nan(size(blockTrialL2R,1),tempTrialIdxMax+tempTrialIdxMin+1)*0.5;

    tempSort = [];
    for i = 1:size(blockTrialL2R,1)
        tempStartIdx = tempTrialIdxMin- (midPoint(i) - blockTrialL2R(i,1));
        tempSort(i) = tempStartIdx;
        tempTrialIdx = blockTrialL2R(i,1):blockTrialL2R(i,4);
        actionTransPlot(i,tempStartIdx+(1:length(tempTrialIdx)) ) = action(tempTrialIdx); 
        stimulusTransPlot(i,tempStartIdx+(1:length(tempTrialIdx)) ) = stimulus(tempTrialIdx); 
        
        actionTransPlotNan(i,tempStartIdx+(1:length(tempTrialIdx)) ) = action(tempTrialIdx); 
        stimulusTransPlotNan(i,tempStartIdx+(1:length(tempTrialIdx)) ) = stimulus(tempTrialIdx); 
    end
    
    [~,sortIdx] = sort(tempSort,'ascend');
    
    plotBin = [-30 30];
    figure; subplot(2,2,1); imagesc(actionTransPlot(sortIdx,:)); 
    xlim([tempTrialIdxMin+plotBin(1) tempTrialIdxMin+plotBin(2)]); caxis([0.2 0.8])
    colormap(redblue); colorbar
    xticks(tempTrialIdxMin+[plotBin(1) 0 plotBin(2)])
    xticklabels(strsplit(int2str([plotBin(1) 0 plotBin(2)])))
    title([mouse 'choice prob']); xlabel('Trials to Transition'); ylabel('nTransitions')
    axis off;
    subplot(2,2,3); hold on;
    ylimm = [0.1 0.9];
    plot([0 0],ylimm,'LineWidth',2,'Color',[0.8 0.8 0.8]);
    plot(-tempTrialIdxMin:tempTrialIdxMax,nanmean(actionTransPlotNan,1),'LineWidth',2,'Color',matlabColors(1));
    xlim(plotBin); ylim(ylimm)

    subplot(2,2,2); imagesc(stimulusTransPlot); 
    xlim([tempTrialIdxMin+plotBin(1) tempTrialIdxMin+plotBin(2)]); caxis([0.2 0.8])
    colormap(redblue); colorbar
    xticks(tempTrialIdxMin+[plotBin(1) 0 plotBin(2)])
    xticklabels(strsplit(int2str([plotBin(1) 0 plotBin(2)])))
    title([mouse 'stimulus prob']); xlabel('Trials to Transition'); ylabel('nTransitions')
    subplot(2,2,4); hold on;
    ylimm = [0.1 0.9];
    plot([0 0],ylimm,'LineWidth',2,'Color',[0.8 0.8 0.8]);
    plot(-tempTrialIdxMin:tempTrialIdxMax,nanmean(stimulusTransPlot,1),'LineWidth',2,'Color',matlabColors(1));
    xlim(plotBin); ylim(ylimm)
    
end

end
