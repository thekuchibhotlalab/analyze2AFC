function plotBiasTransitions(mouseMega)
actionBiasThre = 0.15;
stimCorrThre = 0.3; 
trialBin = 20; corrStim = true;

action = mouseMega.getProp('behav','field','action','matFlag',true);
stimulus = mouseMega.getProp('behav','field','stimulus','matFlag',true);

actL = (action == 1); actL = smoothdata(actL,'movmean',trialBin);
stimL = (stimulus == 1); stimL = smoothdata(stimL,'movmean',trialBin);
[stateFlag, blockL, blockR, normBlock]  = fn_detectBlock(actL,0.5, actionBiasThre,'blockLenThre',1,'blockIntervalThre',5);

if corrStim
    [blockStimFlagL,blockStimProbL] = getBlockStimProb(blockL,stimulus, stimCorrThre);
    [blockStimFlagR,blockStimProbR] = getBlockStimProb(blockR,stimulus, stimCorrThre);
    [blockL, blockR] = correctStimBlock(blockL,blockR,blockStimFlagL,blockStimFlagR);
end


[L2R, R2L] = getTransitions(blockL, blockR);

transTrialThreshold = 30; 

thresholdFlag = (L2R(:,3)-L2R(:,2)) <= transTrialThreshold; L2R = L2R(thresholdFlag,:);
thresholdFlag = (R2L(:,3)-R2L(:,2)) <= transTrialThreshold; R2L = R2L(thresholdFlag,:);

plotTransActionStim(L2R,actL,stimL,[' L2R ']);
plotTransActionStim(R2L,actL,stimL,[' L2R ']);

end
 
function [blockL, blockR] = correctStimBlock(blockL,blockR,blockStimFlagL,blockStimFlagR)
    blockL = fn_structfun(@(x)(x(~blockStimFlagL)),blockL);
    blockR = fn_structfun(@(x)(x(~blockStimFlagR)),blockR);
end

function [A2B, B2A] = getTransitions(idxA, idxB)
    temp = {idxA.start,idxA.end,idxB.start,idxB.end};
    [A2B] = getIdx(idxA.start,idxA.end,idxB.start,idxB.end); % A to B
    [B2A] = getIdx(idxB.start,idxB.end,idxA.start,idxA.end); % B to A

    function [A2B] = getIdx(startA,endA,startB,endB)
        tempDiff = repmat(startB',length(endA),1)-repmat(endA,1,length(startB)); 
        A2B = nan(length(endA),4); 
        for i = 1:length(endA)
            tempIdx = find(tempDiff(i,:) > 0); 
            if ~isempty(tempIdx)
                tempIdx = tempIdx(1);
                A2Bdiff(i) = tempDiff(i,tempIdx);
                A2B(i,1) = startA(i); A2B(i,2) = endA(i); 
                A2B(i,3) = startB(tempIdx); A2B(i,4) = endB(tempIdx);
            end
        end
    end
end

function [blockStimFlag,blockStimProb] = getBlockStimProb(block,stim, probThre)
    nBlock = length(block.start);
    blockStimProb = zeros(size(length(block.start)));
    for i = 1:nBlock 
        blockIdx = block.start(i):block.end(i);
        blockStimProb(i) = mean(stim(blockIdx)==1);
    end
    blockStimFlag = blockStimProb>(0.5+probThre) | blockStimProb<(0.5-probThre); 

end

function plotTransActionStim(blockTrialL2R,action,stimulus,mouse)

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

