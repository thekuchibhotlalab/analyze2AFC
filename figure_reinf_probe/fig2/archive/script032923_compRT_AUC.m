%% LOAD DATA
clear; 
mice ={'zz107','zz109','zz111','zz112','zz113','zz115'};
opsParam.biasBlockType = 'AUC';
allMouse = fn_getObjPT_bin40(mice,opsParam); 
mouseMega = wheel2AFCmega(allMouse);
probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
[label,attributesOrig] = fn_loadWheelCluster();
fastLabel = label'; fastLabel(fastLabel<=2) = 0; fastLabel(fastLabel>2) = 1;
labelCell = {};attributesCell= {};tempTrialSum = [0 cumsum(mouseMega.nTrials)];
for i = 1:mouseMega.nMouse
    labelCell{i} = fastLabel(tempTrialSum(i)+1:tempTrialSum(i+1)); 
    attributesCell{i} = attributesOrig(tempTrialSum(i)+1:tempTrialSum(i+1),:);
    allMouse{i}.behav.label = labelCell{i};

    allMouse{i}.behav.onsetTime = attributesCell{i}(:,1); allMouse{i}.behav.onsetSpeed = attributesCell{i}(:,2); allMouse{i}.behav.onsetDist = attributesCell{i}(:,3); 
    allMouse{i}.behav.totalTime = attributesCell{i}(:,4); allMouse{i}.behav.totalSpeed = attributesCell{i}(:,5); allMouse{i}.behav.totalDist = attributesCell{i}(:,6); 
end
mouseMega = wheel2AFCmega(allMouse);
%[~,badFlag] = getDownSample(mouseMega,[]);

%% 1.0 -- COMPARISON BETWEEN U and LR 
% + is L movement, - is R movement
comp = {'U2L','U2R'};
[allTbl_ctxt] = getTable(allMouse,comp);
%% 1.1 -- raection time comp
makeSixPlot(allTbl_ctxt, 'reactionTime');
makeSixPlot(allTbl_ctxt, 'label');
makeSixPlot(allTbl_ctxt, 'onsetSpeed');
makeSixPlot(allTbl_ctxt, 'onsetTime');
makeSixPlot(allTbl_ctxt, 'onsetDist');
makeSixPlot(allTbl_ctxt, 'totalSpeed');
makeSixPlot(allTbl_ctxt, 'totalTime');
makeSixPlot(allTbl_ctxt, 'totalDist');
makeSixPlot(allTbl_ctxt, 'preMoveAbs');
makeSixPlot(allTbl_ctxt, 'itiMoveAbs');

%%
figure; subplot(1,2,1)
comp = {'U2L','L2U','U2R','R2U'}; a = makeCorrPlot(allMouse,'reactionTime',comp,true);
subplot(1,2,2); comp = {'U2L','L2U','U2R','R2U'}; a = makeCorrPlot(allMouse,'reactionTime',comp,true);
ylim([-0.3 0.3]); xlim([-0.5 0.5])
figure; subplot(1,2,1); comp = {'U2L','L2U','U2R','R2U'}; b = makeCorrPlot(allMouse,'onsetTime',comp,true);
subplot(1,2,2); comp = {'U2L','L2U','U2R','R2U'}; a = makeCorrPlot(allMouse,'onsetTime',comp,true);
ylim([-0.3 0.3]); xlim([-0.5 0.5])
%% total speed and onset time
comp = {'U2L','L2U','U2R','R2U'};
figure; subplot(1,2,1)
comp = {'U2L','L2U','U2R','R2U'}; a = makeCorrPlot(allMouse,'totalSpeed',comp,true);
subplot(1,2,2); comp = {'U2L','L2U','U2R','R2U'}; a = makeCorrPlot(allMouse,'totalSpeed',comp,true);
ylim([-2 2]); xlim([-0.5 0.5])

comp = {'U2L','L2U','U2R','R2U'};
figure; subplot(1,2,1)
comp = {'U2L','L2U','U2R','R2U'}; a = makeCorrPlot(allMouse,'onsetTime',comp,true);
subplot(1,2,2); comp = {'U2L','L2U','U2R','R2U'}; a = makeCorrPlot(allMouse,'onsetTime',comp,true);
ylim([-150 150]); xlim([-0.5 0.5])

%% ITI movement plot
comp = {'U2L','U2R','L2U','R2U'};
[biasDiffL,biasDiffDirL, itemDiffL] = corrPlot(allMouse,'itiMoveL',comp);
[biasDiffR,biasDiffDirR, itemDiffR] = corrPlot(allMouse,'itiMoveR',comp);
[biasDiffLR,biasDiffDirLR, itemDiffLR] = corrPlot(allMouse,'itiMove',comp);
biasDiffDirL = fn_cell2mat(biasDiffDirL,2); itemDiffL = fn_cell2mat(itemDiffL,2);
biasDiffDirR = fn_cell2mat(biasDiffDirR,2); itemDiffR = fn_cell2mat(itemDiffR,2);
biasDiffDirLR = fn_cell2mat(biasDiffDirLR,2); itemDiffLR = fn_cell2mat(itemDiffLR,2);
figure; subplot(1,2,1);
[lm] = fn_plotScatterCorr({biasDiffDirL},{itemDiffL});
xlabel('Increase in leftward bias'); ylabel('Increase of ITI wheel movement in L direction')
subplot(1,2,2); [lm] = fn_plotScatterCorr({biasDiffDirR},{itemDiffR});
xlabel('Increase in leftward bias'); ylabel('Increase of ITI wheel movement in R direction')

figure; subplot(1,2,1);
[lm] = fn_plotScatterCorr({biasDiffDirLR},{itemDiffLR},'scatterColors', {[0.6 0.6 0.6]},'corrPlotArgIn',{'LineWidth',1.5,'Color',[0 0 0]});
xlabel('Increase in left bias'); ylabel('Increase of left ITI movement'); xlim([-1 1]); ylim([-30 30]); yticks(-30:15:30)
subplot(1,2,2);
[lm] = fn_plotScatterCorr({biasDiffDirLR},{itemDiffLR},'scatterColors', {[0.6 0.6 0.6]},'corrPlotArgIn',{'LineWidth',1.5,'Color',[0 0 0]});
xlabel('Increase in left bias'); ylabel('Increase of left ITI movement'); xlim([-0.5 0.5]); ylim([-10 10])


%% label
comp = {'U2L','L2U','U2R','R2U'};
figure; subplot(1,2,1)
comp = {'U2L','L2U','U2R','R2U'}; a = makeCorrPlot(allMouse,'label',comp,true);
subplot(1,2,2); comp = {'U2L','L2U','U2R','R2U'}; a = makeCorrPlot(allMouse,'label',comp,true);
ylim([-2 2]); xlim([-0.5 0.5])

comp = {'U2L','L2U','U2R','R2U'};
figure; subplot(1,2,1)
comp = {'U2L','L2U','U2R','R2U'}; a = makeCorrPlot(allMouse,'onsetTime',comp,true);
subplot(1,2,2); comp = {'U2L','L2U','U2R','R2U'}; a = makeCorrPlot(allMouse,'onsetTime',comp,true);
ylim([-150 150]); xlim([-0.5 0.5])
%% label over learning
tempLabel = {};
for i = 1:length(allMouse); tempLabel{i} = allMouse{i}.behav.label; end
tempLabel = smoothdata(fn_cell2matFillNan(tempLabel),1,'movmean',100,'includenan');

figure; 
[f_line, f_errobar] = fn_plotMeanErrorbar(1:size(tempLabel,1),tempLabel',[0.2 0.2 0.2],[0.6 0.6 0.6],{},{});
%% FUNCTIONS
function itemDiffBalanced = makeCorrPlot(allMouse,itemName,comp,oneFigFlag, colorIn)
    if ~exist('oneFigFlag','var'); oneFigFlag = false; end; if ~exist('colorIn','var'); colorIn = [0.6 0.6 0.6]; end 
    [biasDiff, ~,itemDiff, itemDiffBalanced,itemDiffL,itemDiffR] = corrPlot(allMouse,itemName,comp);

    biasDiff = fn_cell2mat(biasDiff,2); itemDiff = fn_cell2mat(itemDiff,2); itemDiffBalanced = fn_cell2mat(itemDiffBalanced,2);
    itemDiffL = fn_cell2mat(itemDiffL,2); itemDiffR = fn_cell2mat(itemDiffR,2);
    if oneFigFlag
        [lm] = fn_plotScatterCorr({biasDiff},{itemDiffBalanced},'scatterColors', {colorIn}, 'scatterArgIn',{'filled','MarkerFaceAlpha',0.2},'corrPlotArgIn',{'LineWidth',1.5,'Color',[0 0 0]}); 
        title([itemName ', balanced']); xlabel('Change in bias level'); ylabel(['Change in ' itemName]);
    else

        figure; subplot(2,2,1); [lm] = fn_plotScatterCorr({biasDiff},{itemDiff}); title([itemName ', all']);
        xlabel('Change in bias level'); ylabel(['Change in ' itemName])
        subplot(2,2,2); [lm] = fn_plotScatterCorr({biasDiff},{itemDiffBalanced}); title([itemName ', balanced']);
        xlabel('Change in bias level'); ylabel(['Change in ' itemName])
        subplot(2,2,3); [lm] = fn_plotScatterCorr({biasDiff},{itemDiffL}); title([itemName ', L']);
        xlabel('Change in bias level'); ylabel(['Change in ' itemName])
        subplot(2,2,4); [lm] = fn_plotScatterCorr({biasDiff},{itemDiffR}); title([itemName ', R']);
        xlabel('Change in bias level'); ylabel(['Change in ' itemName])
    end 
end


function makeSixPlot(allTbl_ctxt, itemName)
figure; 
subplot(2,3,1);makePlot(allTbl_ctxt, itemName, 'balanceLR', 'fast'); title([itemName ' balanceLR -- Fast'])
xticks([1 2]); xticklabels({'unbiased','biased'})
subplot(2,3,2);makePlot(allTbl_ctxt, itemName, 'balanceLR', 'slow'); title([itemName ' balanceLR -- slow'])
xticks([1 2]); xticklabels({'unbiased','biased'})
subplot(2,3,3);makePlot(allTbl_ctxt, itemName, 'balanceLR', 'all'); title([itemName ' balanceLR -- all'])
xticks([1 2]); xticklabels({'unbiased','biased'})
subplot(2,3,4);makePlot(allTbl_ctxt, itemName, 'allTrial', 'fast'); title([itemName ' allTrial -- Fast'])
xticks([1 2]); xticklabels({'unbiased','biased'})
subplot(2,3,5);makePlot(allTbl_ctxt, itemName, 'allTrial', 'slow'); title([itemName ' allTrial -- slow'])
xticks([1 2]); xticklabels({'unbiased','biased'})
subplot(2,3,6);makePlot(allTbl_ctxt, itemName, 'allTrial', 'all'); title([itemName ' allTrial -- all'])
xticks([1 2]); xticklabels({'unbiased','biased'})
end

function makePlot(allTbl_ctxt, itemName, compType, label)
    switch compType
    case 'sepLR'
        [item] = getTableItemLR(allTbl_ctxt,itemName,label); % label is all or fast or slow
        p = fn_plotComparison({item.unbias(:,1),item.bias(:,1)},'paired',true,'compType','errorbarWithDot');

    case 'allTrial'
        [item] = getTableItemAllTrials(allTbl_ctxt,itemName,label); % label is all or fast or slow
        p = fn_plotComparison({item.unbias(:,1),item.bias(:,1)},'paired',true,'compType','errorbarWithDot');
    case 'balanceLR'
        [item] = getTableItem(allTbl_ctxt,itemName,label); % label is all or fast or slow
        p = fn_plotComparison({item.unbias(:,1),item.bias(:,1)},'paired',true,'compType','errorbarWithDot');
    end 



end

function [item] = getTableItem(allTbl_ctxt,itemName,label)
item.unbias = zeros(length(allTbl_ctxt),1);item.bias = zeros(length(allTbl_ctxt),1);
for i = 1:length(allTbl_ctxt)
    tempUnbias = nan(2,1); tempBias = nan(2,1);
    for dir = 1:2 % loop through action type      
        if strcmp(label,'fast'); tempIdx =  allTbl_ctxt{i}.action == dir & allTbl_ctxt{i}.label == 0;
        elseif strcmp(label,'slow'); tempIdx =  allTbl_ctxt{i}.action == dir & allTbl_ctxt{i}.label == 1;
        else tempIdx =  allTbl_ctxt{i}.action == dir; end 
        tempUnbias(dir) =  nanmean(allTbl_ctxt{i}.(itemName)(tempIdx & allTbl_ctxt{i}.biasLabel==0));
        tempBias(dir) =  nanmean(allTbl_ctxt{i}.(itemName)(tempIdx & allTbl_ctxt{i}.biasLabel==1));        
    end
    item.unbias(i)= nanmean(tempUnbias); item.bias(i) = nanmean(tempBias); 
end
end

function [item] = getTableItemLR(allTbl_ctxt,itemName,label)
item.unbias = zeros(length(allTbl_ctxt),2); item.bias = zeros(length(allTbl_ctxt),2);
for i = 1:length(allTbl_ctxt)
    tempUnbias = nan(2,1); tempBias= nan(2,1);
    for dir = 1:2 % loop through action type
        if strcmp(label,'fast'); tempIdx =  allTbl_ctxt{i}.action == dir & allTbl_ctxt{i}.label == 0;
        elseif strcmp(label,'slow'); tempIdx =  allTbl_ctxt{i}.action == dir & allTbl_ctxt{i}.label == 1;
        else tempIdx =  allTbl_ctxt{i}.action == dir; end 
        tempUnbias(dir) =  nanmean(allTbl_ctxt{i}.(itemName)(tempIdx & allTbl_ctxt{i}.biasLabel==0));
        tempBias(dir) =  nanmean(allTbl_ctxt{i}.(itemName)(tempIdx & allTbl_ctxt{i}.biasLabel==1));          
    end
    item.unbias(i,:)= tempUnbias; item.bias(i,:) = tempBias; 
end
end

function [item] = getTableItemAllTrials(allTbl_ctxt,itemName,label)
item.unbias = zeros(length(allTbl_ctxt),1);item.bias = zeros(length(allTbl_ctxt),1);
for i = 1:length(allTbl_ctxt)
    if strcmp(label,'fast'); tempIdx =  allTbl_ctxt{i}.label == 0;
    elseif strcmp(label,'slow'); tempIdx =  allTbl_ctxt{i}.label == 1;
    else tempIdx = true(size(allTbl_ctxt{i}.label)); end 
    tempUnbias = nanmean(allTbl_ctxt{i}.(itemName)(tempIdx&allTbl_ctxt{i}.biasLabel==0));
    tempBias = nanmean(allTbl_ctxt{i}.(itemName)(tempIdx&allTbl_ctxt{i}.biasLabel==1));
    item.unbias(i,:)= tempUnbias; item.bias(i,:) = tempBias; 
end
end

function [allTbl_ctxt] = getTable(allMouse,comp)

nMouse = length(allMouse); 
allTbl_ctxt = cell(1,nMouse); countL = cell(1,nMouse); countR = cell(1,nMouse); countBlock = cell(1,nMouse);
for i = 1:nMouse
    allTbl_ctxt{i} = allMouse{i}.behav; 
    accThre = [0.7 nan]; windowSize= 400; 
    [~, learningOnset] = allMouse{i}.fn_findLearningOnset(accThre,windowSize);
    selTrial = learningOnset-300:learningOnset+300;
    nTrials = size(allMouse{i}.behav,1); 
    selTrial(selTrial<=0) = []; selTrial(selTrial>nTrials) = [];
    %dayStart = allMouse{i}.behav.day(selTrial(1)); dayEnd = allMouse{i}.behav.day(selTrial(end));
    befAft = 10;
    blockStartFlag = contains(allMouse{i}.biasBlock.transID,comp);
    blockStart = (allMouse{i}.biasBlock.trans(blockStartFlag,2)); blockEnd = (allMouse{i}.biasBlock.trans(blockStartFlag,4));
    disp(length(blockStart))
    blockBefIdx = cell(1,length(blockStart)); blockAftIdx = cell(1,length(blockStart));
   
    for j = 1:length(blockStart)
        tempBefCount = min(blockStart(j)-1,befAft); tempAftCount = min(nTrials-blockEnd(j),befAft);
        tempCount = min([tempBefCount, tempAftCount]); 
        [maxBias,maxBiasIdx] = max(abs(allMouse{i}.behav.bias(blockStart(j):blockEnd(j))));
        tempIdx = blockStart(j) + maxBiasIdx - befAft + 1 : blockStart(j) + maxBiasIdx + befAft;        
        if maxBias >= 0.4
            blockBefIdx{j} = [blockStart(j)-tempCount:blockStart(j)-1 blockEnd(j)+1:blockEnd(j)+tempCount];
            %blockBefIdx{j} = [blockStart(j)+1:tempCount+blockStart(j) blockEnd(j)-tempCount+1:blockEnd(j)];
            blockAftIdx{j} = tempIdx;
        else 
            blockBefIdx{j} = [];
            blockAftIdx{j} = [];
        end
        countL{i}(j,:) = [sum(allTbl_ctxt{i}.action(blockBefIdx{j})==1) sum(allTbl_ctxt{i}.action(blockAftIdx{j})==1)]; 
        countBlock{i}(j) = blockStart(j);
        countR{i}(j,:) = [sum(allTbl_ctxt{i}.action(blockBefIdx{j})==2) sum(allTbl_ctxt{i}.action(blockAftIdx{j})==2)]; 

    end 
    blockBefIdx = fn_cell2mat(blockBefIdx,2);
    blockAftIdx = fn_cell2mat(blockAftIdx,2);

    tempBiasLabel = nan(size(allMouse{i}.behav,1),1); tempBiasLabel(blockBefIdx) = 0; tempBiasLabel(blockAftIdx) = 1;
    allTbl_ctxt{i}.biasLabel = tempBiasLabel;
    tempSelIdx = [blockBefIdx(:);blockAftIdx(:)];
    allTbl_ctxt{i} = allTbl_ctxt{i}(tempSelIdx,:);
end 
end

function  [biasDiff,biasDiffDir, itemDiff,itemDiffBalanced,itemDiffL,itemDiffR] = corrPlot(allMouse,item,comp)
nMouse = length(allMouse); 
allTbl_ctxt = cell(1,nMouse); countL = cell(1,nMouse); countR = cell(1,nMouse); countBlock = cell(1,nMouse);
for i = 1:nMouse
    allTbl_ctxt{i} = allMouse{i}.behav; 
    nTrials = size(allMouse{i}.behav,1); 

    binSize = 40; moveBin = 20;  
    blockStartFlag = contains(allMouse{i}.biasBlock.transID,comp);
    blockStart = (allMouse{i}.biasBlock.trans(blockStartFlag,2)); 
    blockEnd = (allMouse{i}.biasBlock.trans(blockStartFlag,4));

    biasDiffDir{i} = []; biasDiff{i} = []; itemDiff{i} = []; itemDiffBalanced{i} = []; itemDiffL{i} = [];itemDiffR{i} = [];
    for j =1:length(blockStart)
        tempLen = blockEnd(j) - blockStart(j)+1;
        if tempLen >= binSize*2
            tempCrop = ceil(mod(tempLen,moveBin)/2); 
            tempStart = blockStart(j)+tempCrop; tempEnd = tempStart+floor(tempLen/moveBin)*moveBin-1;
            tempBin = tempStart:moveBin:tempEnd; tempMaxIdx = find(tempBin == (max(tempBin) - binSize*2));
            for k = 1:tempMaxIdx
                idx1 = tempBin(k):tempBin(k)+binSize-1; idx2 = tempBin(k)+binSize:tempBin(k)+binSize*2-1;
                [tempBias1,~,~,~,~,~] = fn_getAccBias(allMouse{i}.behav.stimulus(idx1),...
                    allMouse{i}.behav.responseType(idx1)==1);
                [tempBias2,~,~,~,~,~] = fn_getAccBias(allMouse{i}.behav.stimulus(idx2),...
                    allMouse{i}.behav.responseType(idx2)==1);
                tempBiasDiff = abs(tempBias2) - abs(tempBias1);
                tempBiasDiffDir = (tempBias2) - (tempBias1);
                tempItemDiff = nanmean(allMouse{i}.behav.(item)(idx2))-nanmean(allMouse{i}.behav.(item)(idx1));
                
                tempTable1 = allMouse{i}.behav(idx1,:); tempTable2 = allMouse{i}.behav(idx2,:);
                tempItem1 = nanmean(tempTable1.(item)(tempTable1.action==1))/2 + nanmean(tempTable1.(item)(tempTable1.action==2))/2;
                tempItem2 = nanmean(tempTable2.(item)(tempTable2.action==1))/2 + nanmean(tempTable2.(item)(tempTable2.action==2))/2;
                tempItemDiffBalanced = tempItem2 - tempItem1;
                tempItemDiffL = nanmean(tempTable2.(item)(tempTable2.action==1)) - nanmean(tempTable1.(item)(tempTable1.action==1));
                tempItemDiffR = nanmean(tempTable2.(item)(tempTable2.action==2)) - nanmean(tempTable1.(item)(tempTable1.action==2));
                biasDiff{i} = [biasDiff{i} tempBiasDiff];  biasDiffDir{i} = [biasDiffDir{i} tempBiasDiffDir]; itemDiff{i} = [itemDiff{i} tempItemDiff]; 
                itemDiffL{i} = [itemDiffL{i} tempItemDiffL];itemDiffR{i} = [itemDiffR{i} tempItemDiffR];
                itemDiffBalanced{i} = [itemDiffBalanced{i} tempItemDiffBalanced];
            end

        else

        end 
    end

end

end