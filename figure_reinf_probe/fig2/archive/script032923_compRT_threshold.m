%% LOAD DATA
clear; 
mice ={'zz107','zz109','zz111','zz112','zz113','zz115'};
opsParam.biasBlockType = 'threshold';
allMouse = fn_getObjPT_bin30(mice,opsParam); 
mouseMega = wheel2AFCmega(allMouse);
probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
[label,attributesOrig] = fn_loadWheelCluster();
fastLabel = label'; fastLabel(fastLabel<=2) = 0; fastLabel(fastLabel>2) = 1;
%fastLabel(fastLabel==2) = 0;fastLabel(fastLabel==3) = 0;
%fastLabel(fastLabel==1) = 1; fastLabel(fastLabel==4) = 1;
labelCell = {};attributesCell= {};tempTrialSum = [0 cumsum(mouseMega.nTrials)];
for i = 1:mouseMega.nMouse
    labelCell{i} = fastLabel(tempTrialSum(i)+1:tempTrialSum(i+1)); 
    attributesCell{i} = attributesOrig(tempTrialSum(i)+1:tempTrialSum(i+1),:);
    allMouse{i}.behav.label = labelCell{i};
end
mouseMega = wheel2AFCmega(allMouse);
%[~,badFlag] = getDownSample(mouseMega,[]);

%% 1.0 -- COMPARISON BETWEEN U and LR 
% + is L movement, - is R movement
comp = {'U2L','U2R'};
%comp = {'L2U','R2U'};
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

%% 2.0 -- COMPARISON BETWEEN U and L
% + is L movement, - is R movement
comp = {'U2R'};
%comp = {'L2U','R2U'};
[allTbl_ctxt] = getTableItemLR(allMouse,comp);
%% 2.1 -- COMPARISON BETWEEN U and L

makeSixPlot(allTbl_ctxt, 'itiMoveL');
%%
[item] = getTableItemAllTrials(allTbl_ctxt,'preMove',false);
%%
[a] = getTableItem(allTbl_ctxt,'totalDist',false);
figure; p = fn_plotComparison({a.unbias(:,1),a.bias(:,1)},'paired',true,'compType','errorbarWithDot');

%% NEW CODE TO CONFIRM COMPARISON

[item] = getTableItemLR(allTbl_ctxt,'preMoveL','all');
%% FUNCTIONS
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

    befAft = 20;
    blockStartFlag = contains(allMouse{i}.biasBlock.transID,comp);
    blockStart = (allMouse{i}.biasBlock.trans(blockStartFlag,2));
    disp(length(blockStart))
    blockBefIdx = cell(1,length(blockStart)); blockAftIdx = cell(1,length(blockStart));

    
    for j = 1:length(blockStart)
        tempBefCount = min(blockStart(j)-1,befAft); tempAftCount = min(nTrials-blockStart(j)+1,befAft);
        tempCount = min([tempBefCount, tempAftCount]);
        blockBefIdx{j} = blockStart(j)-tempCount:blockStart(j)-1;
        blockAftIdx{j} = blockStart(j):blockStart(j)+tempCount-1;

        
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