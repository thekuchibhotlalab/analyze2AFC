%% LOAD DATA
clear; 
mice ={'zz107','zz109','zz111','zz112','zz113','zz115'};
opsParam.biasBlockType = 'threshold';
allMouse = fn_getObjPT_bin30(mice,opsParam); mouseMega = wheel2AFCmega(allMouse);

probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
[label,attributesOrig] = fn_loadWheelCluster();
fastLabel = label'; fastLabel(fastLabel<=2) = 0; fastLabel(fastLabel>2) = 1;
labelCell = {};attributesCell= {};tempTrialSum = [0 cumsum(mouseMega.nTrials)];
for i = 1:mouseMega.nMouse
    labelCell{i} = fastLabel(tempTrialSum(i)+1:tempTrialSum(i+1)); 
    attributesCell{i} = attributesOrig(tempTrialSum(i)+1:tempTrialSum(i+1),:); 
end
[~,badFlag] = getDownSample(mouseMega,[]);

%% NEW CODE TO CONFIRM COMPARISON
allTbl_ctxt = {}; wheelAttrName = {'startTime','startSpeed','startDist','totalTime','totalSpeed','totalDist'};
biasRT_all_by_block = cell(1,length(mouseMega.nMouse));
for i = 1:mouseMega.nMouse
    accThre = [0.7 nan]; windowSize= 400; 
    [~, learningOnset] = allMouse{i}.fn_findLearningOnset(accThre,windowSize);
    selTrial = learningOnset-300:learningOnset+300;
    nTrials = size(allMouse{i}.behav,1); 
    selTrial(selTrial<=0) = []; selTrial(selTrial>nTrials) = [];
    dayStart = allMouse{i}.behav.day(selTrial(1)); dayEnd = allMouse{i}.behav.day(selTrial(end));
    %startTrial =  find(allMouse{i}.behav.day==dayStart,1); endTrial =  find(allMouse{i}.behav.day==dayEnd,1,'last');
    %allBlockStart = [allMouse{i}.biasBlock.blockL.start;allMouse{i}.biasBlock.blockR.start];

    tempSpeed = diff(allMouse{i}.behavVar.wheelPos_aligned,1,2);
    preMoveAbs = nansum(abs(tempSpeed(:,1:2000)),2);
    preMove = nansum((tempSpeed(:,1:2000)),2);
    allTbl_ctxt{i} = allMouse{i}.behav;
    allTbl_ctxt{i}.label = labelCell{i};
    %allTbl_ctxt{i}.biasLabel = tempBiasLabel(tempSelIdx);
    for j = 1:length(wheelAttrName); allTbl_ctxt{i}.(wheelAttrName{j}) = attributesCell{i}(:,j); end

    allTbl_ctxt{i}.preMoveAbs = preMoveAbs;
    allTbl_ctxt{i}.preMove = preMove .* (allTbl_ctxt{i}.action*2-3);

    allTbl_ctxt{i}.reactionTime(allTbl_ctxt{i}.reactionTime>=2.5) = nan; 
    allTbl_ctxt{i}.reactionTime(allTbl_ctxt{i}.reactionTime<=0) = nan; 
    

    befAft = 20;
    blockStartFlag = contains(allMouse{i}.biasBlock.transID,{'U2L','U2R'});
    blockStart = (allMouse{i}.biasBlock.trans(blockStartFlag,2));
    disp(length(blockStart))
    blockBefIdx = cell(1,length(blockStart)); blockAftIdx = cell(1,length(blockStart));
    for j = 1:length(blockStart)
        tempBefCount = min(blockStart(j)-1,befAft); tempAftCount = min(nTrials-blockStart(j)+1,befAft);
        tempCount = min([tempBefCount, tempAftCount]);
        blockBefIdx{j} = blockStart(j)-tempCount:blockStart(j)-1;
        blockAftIdx{j} = blockStart(j):blockStart(j)+tempCount-1;
        tempBiasLabel = nan(size(allMouse{i}.behav,1),1); tempBiasLabel(blockBefIdx{j}) = 0; tempBiasLabel(blockAftIdx{j}) = 1;
        tempTbl = allTbl_ctxt{i}; tempTbl.biasLabel = tempBiasLabel;
        [tempBiasRT_all] = getTableItem({tempTbl},'reactionTime',false);
        biasRT_all_by_block{i}(j,:) =  [tempBiasRT_all.unbias tempBiasRT_all.bias];

        [tempBiasRT_all] = getTableItem({tempTbl},'reactionTime',false);
        biasRT_all_by_block{i}(j,:) =  [tempBiasRT_all.unbias tempBiasRT_all.bias];
    end 
    blockBefIdx = fn_cell2mat(blockBefIdx,2);
    blockAftIdx = fn_cell2mat(blockAftIdx,2);
    tempBiasLabel = nan(size(allMouse{i}.behav,1),1); tempBiasLabel(blockBefIdx) = 0; tempBiasLabel(blockAftIdx) = 1;
    tempSelIdx = [blockBefIdx(:);blockAftIdx(:)];   
end

[biasRT] = getTableItem(allTbl_ctxt,'reactionTime',true);
[biasRT_all] = getTableItem(allTbl_ctxt,'reactionTime',false);
[biasLabel] = getTableItem(allTbl_ctxt,'label',false);
[biasPreMove] = getTableItem(allTbl_ctxt,'preMove',false);
biasWheel = {}; 
for i = 1:length(wheelAttrName); biasWheel{i} = getTableItem(allTbl_ctxt,wheelAttrName{i},false); end

%%
figure; subplot(2,4,1)
p = fn_plotComparison({biasRT.unbias(:,1),biasRT.bias(:,1)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks(1:2); ylim([0.3 0.7]); yticks(0.3:0.2:0.7);title(['p=' num2str(p,'%.4f')])
subplot(2,4,2)
p = fn_plotComparison({biasRT_all.unbias(:,1),biasRT_all.bias(:,1)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks(1:2); ylim([0.3 1.1]); yticks(0.3:0.4:1.1);title(['p=' num2str(p,'%.4f')])

ylimCell = {[0 800],[3 6],[0.8 1.2],[0 500],[3 10],[1 1.5]};
for i = 1:length(biasWheel)
    subplot(2,4,i+2);
    p = fn_plotComparison({biasWheel{i}.unbias,biasWheel{i}.bias},'paired',true,'compType','errorbarWithDot');
    title(['p=' num2str(p,'%.4f')])
    xlim([0.6 2.4]); xticks(1:2); ylim(ylimCell{i}); 
end
figure;
p = fn_plotComparison({biasPreMove.unbias,biasPreMove.bias},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks(1:3);title(['p=' num2str(p,'%.4f')])
%% session by session 
blockComp_RT = fn_cell2mat(cellfun(@(x)(mean(x,1)),biasRT_all_by_block,'UniformOutput',false),1);
figure; p=fn_plotComparison({blockComp_RT(:,1),blockComp_RT(:,2)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks(1:2);title(['p=' num2str(p,'%.4f')])
%% FUNCTIONS
function [item] = getTableItem(allTbl_ctxt,itemName,separateLabel)
if nargin == 2; separateLabel = false; end
if separateLabel
    item.unbias = zeros(length(allTbl_ctxt),2); item.bias = zeros(length(allTbl_ctxt),2);
else
    item.unbias = zeros(length(allTbl_ctxt),1);item.bias = zeros(length(allTbl_ctxt),1);
end
for i = 1:length(allTbl_ctxt)
    if separateLabel
        tempUnbias = nan(2,2); tempBias= nan(2,2);
    else
        tempUnbias = nan(2,1); tempBias = nan(2,1);
    end

    for dir = 1:2 % loop through action type      
        if separateLabel
            for l = 1:2 % loop through label
                tempIdx =  allTbl_ctxt{i}.action == dir & allTbl_ctxt{i}.label == (l-1);
                tempUnbias(dir,l) =  nanmean(allTbl_ctxt{i}.(itemName)(tempIdx & allTbl_ctxt{i}.biasLabel==0));
                tempBias(dir,l) =  nanmean(allTbl_ctxt{i}.(itemName)(tempIdx & allTbl_ctxt{i}.biasLabel==1));
            end
        else
            tempUnbias(dir) = mean(allTbl_ctxt{i}.(itemName)(allTbl_ctxt{i}.action == dir & allTbl_ctxt{i}.biasLabel==0));
            tempBias(dir) = mean(allTbl_ctxt{i}.(itemName)(allTbl_ctxt{i}.action == dir & allTbl_ctxt{i}.biasLabel==1));
        end  
    end
    item.unbias(i,:)= nanmean(tempUnbias,1); item.bias(i,:) = nanmean(tempBias,1); 
end

end
