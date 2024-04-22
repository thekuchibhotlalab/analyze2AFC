%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
%mice ={'zz107','zz109','zz111','zz112','zz113','zz115'};

opsParam.biasBlockType = 'threshold';
allMouse = fn_getObjPT_bin30(mice,opsParam); mouseMega = wheel2AFCmega(allMouse);

probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;

%% NEW CODE TO CONFIRM COMPARISON
allTbl_ctxt = {}; aR_allani = {mouseMega.nMouse,2};
for i = 1:mouseMega.nMouse
    accThre = [0.7 nan]; windowSize= 400; 
    [~, learningOnset] = allMouse{i}.fn_findLearningOnset(accThre,windowSize);
    selTrial = learningOnset-300:learningOnset+300;
    nTrials = size(allMouse{i}.behav,1); 
    selTrial(selTrial<=0) = []; selTrial(selTrial>nTrials) = [];
    dayStart = allMouse{i}.behav.day(selTrial(1)); dayEnd = allMouse{i}.behav.day(selTrial(end));

    befAft = 30;

    blockStartFlag = contains(allMouse{i}.biasBlock.transID,{'U2L','U2R'});
    blockStart = (allMouse{i}.biasBlock.trans(blockStartFlag,2));

    blockBefIdx = cell(1,length(blockStart)); blockAftIdx = cell(1,length(blockStart));
    for j = 1:length(blockStart)
        tempBefCount = min(blockStart(j)-1,befAft); tempAftCount = min(nTrials-blockStart(j)+1,befAft);
        tempCount = min([tempBefCount, tempAftCount]);
        blockBefIdx{j} = blockStart(j)-tempCount:blockStart(j)-1;
        blockAftIdx{j} = blockStart(j):blockStart(j)+tempCount-1;
        % GET ACTION RATE DATA FOR EACH BLOCK
        tempMasterTrial = allMouse{i}.behav.masterTrialnum(blockBefIdx{j});
        tempTrial = allMouse{i}.behav.trialnum(blockBefIdx{j});
        aR_allani{i,1}(j) = (tempTrial(end) - tempTrial(1)+1) / (tempMasterTrial(end) - tempMasterTrial(1)+1);

        tempMasterTrial = allMouse{i}.behav.masterTrialnum(blockAftIdx{j});
        tempTrial = allMouse{i}.behav.trialnum(blockAftIdx{j});
        aR_allani{i,2}(j) = (tempTrial(end) - tempTrial(1)+1) ./ (tempMasterTrial(end) - tempMasterTrial(1)+1);
    end 

    blockBefIdx = fn_cell2mat(blockBefIdx,2);
    blockAftIdx = fn_cell2mat(blockAftIdx,2);

    tempBiasLabel = nan(size(allMouse{i}.behav,1),1); tempBiasLabel(blockBefIdx) = 0; tempBiasLabel(blockAftIdx) = 1;
    tempSelIdx = [blockBefIdx(:);blockAftIdx(:)];
    tempSelIdx(tempSelIdx>size(allMouse{i}.behav,1)) = [];

    allTbl_ctxt{i} = allMouse{i}.behav(tempSelIdx,:);
    allTbl_ctxt{i}.biasLabel = tempBiasLabel(tempSelIdx);
    allTbl_ctxt{i}.reactionTime(allTbl_ctxt{i}.reactionTime>=2.5) = nan; 
    allTbl_ctxt{i}.reactionTime(allTbl_ctxt{i}.reactionTime<=0) = nan; 
    
end

[biasRT] = getTableItem(allTbl_ctxt,'reactionTime');

%%
figure; 
p = fn_plotComparison({biasRT.unbias(:,1),biasRT.bias(:,1)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks(1:2); ylim([0.3 1.2]); yticks(0.3:0.3:1.2);title(['p=' num2str(p,'%.4f')])

%% Look at action rate
figure;
aR = cellfun(@(x)(mean(x)),aR_allani,'UniformOutput',true);
p = fn_plotComparison({aR(:,1),aR(:,2)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks(1:2); ylim([0.3 1.2]); yticks(0.3:0.3:1.2);title(['p=' num2str(p,'%.4f')])
xticklabels({'unbias','bias'})


%% FUNCTIONS
function [item] = getTableItem(allTbl_ctxt,itemName)

item.unbias = zeros(length(allTbl_ctxt),1);item.bias = zeros(length(allTbl_ctxt),1);

for i = 1:length(allTbl_ctxt)
    tempUnbias = nan(2,1); tempBias = nan(2,1);
    for dir = 1:2 % loop through action type      
        tempUnbias(dir) = mean(allTbl_ctxt{i}.(itemName)(allTbl_ctxt{i}.action == dir & allTbl_ctxt{i}.biasLabel==0));
        tempBias(dir) = mean(allTbl_ctxt{i}.(itemName)(allTbl_ctxt{i}.action == dir & allTbl_ctxt{i}.biasLabel==1));
    end
    item.unbias(i,:)= nanmean(tempUnbias,1); item.bias(i,:) = nanmean(tempBias,1); 
end
end
