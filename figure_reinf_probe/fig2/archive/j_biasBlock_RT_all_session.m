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
countL = cell(1,length(mouseMega.nMouse)); countR = cell(1,length(mouseMega.nMouse)); countBlock = cell(1,length(mouseMega.nMouse));
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
    allTbl_ctxt{i} = allMouse{i}.behav;  allTbl_ctxt{i}.trialNum = (1:size(allMouse{i}.behav,1))';
    allTbl_ctxt{i}.label = labelCell{i};
    %allTbl_ctxt{i}.biasLabel = tempBiasLabel(tempSelIdx);
    for j = 1:length(wheelAttrName); allTbl_ctxt{i}.(wheelAttrName{j}) = attributesCell{i}(:,j); end

    allTbl_ctxt{i}.preMoveAbs = preMoveAbs;
    allTbl_ctxt{i}.preMove = preMove .* (allTbl_ctxt{i}.action*2-3);

    allTbl_ctxt{i}.reactionTime(allTbl_ctxt{i}.reactionTime>=2.5) = nan; 
    allTbl_ctxt{i}.reactionTime(allTbl_ctxt{i}.reactionTime<=0) = nan; 
    

    befAft = 20;
    blockStartFlag = contains(allMouse{i}.biasBlock.transID,{'U2L'});
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



%%
[biasRT] = getTableItem(allTbl_ctxt,'reactionTime',true);
[biasRT_all] = getTableItem(allTbl_ctxt,'reactionTime',false);
[biasLabel] = getTableItem(allTbl_ctxt,'label',false);
[biasPreMove] = getTableItem(allTbl_ctxt,'preMove',false);
biasWheel = {}; 
for i = 1:length(wheelAttrName); biasWheel{i} = getTableItem(allTbl_ctxt,wheelAttrName{i},false); end
figure; subplot(2,4,1)
p = fn_plotComparison({biasRT.unbias(:,1),biasRT.bias(:,1)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks(1:2); ylim([0.3 0.7]); yticks(0.3:0.2:0.7);title(['p=' num2str(p,'%.4f')]);xticklabels({'unbias','bias'});
subplot(2,4,2)
p = fn_plotComparison({biasRT_all.unbias(:,1),biasRT_all.bias(:,1)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks(1:2); ylim([0.3 1.1]); yticks(0.3:0.4:1.1);title(['p=' num2str(p,'%.4f')]);xticklabels({'unbias','bias'});

ylimCell = {[0 800],[3 6],[0.8 1.2],[0 500],[3 10],[1 1.5]};
for i = 1:length(biasWheel)
    subplot(2,4,i+2);
    p = fn_plotComparison({biasWheel{i}.unbias,biasWheel{i}.bias},'paired',true,'compType','errorbarWithDot');
    title(['p=' num2str(p,'%.4f')])
    xlim([0.6 2.4]); xticks(1:2); xticklabels({'bias','unbias'});ylim(ylimCell{i}); 
end
figure;
p = fn_plotComparison({biasLabel.unbias,biasLabel.bias},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks(1:3);title(['p=' num2str(p,'%.4f')]); xticklabels({'unbias','bias'});

figure;
p = fn_plotComparison({biasPreMove.unbias,biasPreMove.bias},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks(1:3);title(['p=' num2str(p,'%.4f')]); xticklabels({'unbias','bias'});
%% snumber of trials
figure; 

countBlockmat = fn_cell2mat(countBlock,2);
Rmat = fn_cell2mat(countR,1); 
Lmat = fn_cell2mat(countL,1);
[binCount,b,c] = histcounts(countBlockmat);

nBins = length(binCount);
Rmat_count = zeros(nBins,2); Lmat_count = zeros(nBins,2);
for i = 1:length(c)
    Lmat_count(c(i),:) = Lmat_count(c(i),:) + Lmat(i,:);
    Rmat_count(c(i),:) = Rmat_count(c(i),:) + Rmat(i,:);
end
subplot(1,2,1); plot(b(2:end)-250,Lmat_count); legend({'unbiased','biased'}); 
xlabel('Trials in training'); ylabel('trials taken for comparison'); title('Left')
subplot(1,2,2); plot(b(2:end)-250,Rmat_count); legend({'unbiased','biased'}); 
xlabel('Trials in training'); ylabel('trials taken for comparison'); title('Right')

%%
befMotor = cell(mouseMega.nMouse,2); aftMotor = cell(mouseMega.nMouse,2);
befWheel = cell(mouseMega.nMouse,2); aftWheel = cell(mouseMega.nMouse,2);
for i = 1:mouseMega.nMouse
    befTrialsL = allTbl_ctxt{i}.trialNum(allTbl_ctxt{i}.biasLabel==0 & allTbl_ctxt{i}.action==1);
    befTrialsR = allTbl_ctxt{i}.trialNum(allTbl_ctxt{i}.biasLabel==0 & allTbl_ctxt{i}.action==2);

    aftTrialsL = allTbl_ctxt{i}.trialNum(allTbl_ctxt{i}.biasLabel==1 & allTbl_ctxt{i}.action==1);
    aftTrialsR = allTbl_ctxt{i}.trialNum(allTbl_ctxt{i}.biasLabel==1 & allTbl_ctxt{i}.action==2);
    befMotor{i,1} = allMouse{i}.behavVar.wheelPos_aligned(befTrialsL,:);
    befMotor{i,2} = allMouse{i}.behavVar.wheelPos_aligned(befTrialsR,:);
    aftMotor{i,1} = allMouse{i}.behavVar.wheelPos_aligned(aftTrialsL,:);
    aftMotor{i,2} = allMouse{i}.behavVar.wheelPos_aligned(aftTrialsR,:);

    befWheel{i,1} = allTbl_ctxt{i}{allTbl_ctxt{i}.biasLabel==0 & allTbl_ctxt{i}.action==1,42:47};
    befWheel{i,2} = allTbl_ctxt{i}{allTbl_ctxt{i}.biasLabel==0 & allTbl_ctxt{i}.action==2,42:47};
    aftWheel{i,1} = allTbl_ctxt{i}{allTbl_ctxt{i}.biasLabel==1 & allTbl_ctxt{i}.action==1,42:47};
    aftWheel{i,2} = allTbl_ctxt{i}{allTbl_ctxt{i}.biasLabel==1 & allTbl_ctxt{i}.action==2,42:47};
end
befSpeedAvg = cellfun(@(x)(mean(diff((x),1,2),1)),befMotor,'UniformOutput',false);
aftSpeedAvg = cellfun(@(x)(mean(diff((x),1,2),1)),aftMotor,'UniformOutput',false);

befMotorAvg = cellfun(@(x)(mean((x),1)),befMotor,'UniformOutput',false);
aftMotorAvg = cellfun(@(x)(mean((x),1)),aftMotor,'UniformOutput',false);

figure; 
for i = 1:mouseMega.nMouse
    subplot(2,3,i); 
    plot(-3:0.001:5,mean(befMotorAvg{i,1},1)); hold on; plot(-3:0.001:5,mean(aftMotorAvg{i,1},1)); 
    xlabel('Time from sound onset'); ylabel('Position toward left'); title(['mouse' int2str(i)])
end
figure; 
for i = 1:mouseMega.nMouse
    subplot(2,3,i); 
    plot(-3:0.001:5,mean(befMotorAvg{i,2},1)); hold on; plot(-3:0.001:5,mean(aftMotorAvg{i,2},1)); 
    xlabel('Time from sound onset'); ylabel('Position toward left'); title(['mouse' int2str(i)])
end

smoothWindow = 50;
figure; 
for i = 1:mouseMega.nMouse
    subplot(2,3,i); 
    plot(-3+0.001:0.001:5,smoothdata(befSpeedAvg{i,1},'movmean',smoothWindow)); hold on; plot(-3+0.001:0.001:5,smoothdata(aftSpeedAvg{i,1},'movmean',smoothWindow)); 
    xlabel('Time from sound onset'); ylabel('Speed toward left'); title(['mouse' int2str(i)])
end
figure; 
for i = 1:mouseMega.nMouse
    subplot(2,3,i); 
    plot(-3+0.001:0.001:5,smoothdata(befSpeedAvg{i,2},'movmean',smoothWindow)); hold on; plot(-3+0.001:0.001:5,smoothdata(aftSpeedAvg{i,2},'movmean',smoothWindow)); 
    xlabel('Time from sound onset'); ylabel('Speed toward left'); title(['mouse' int2str(i)])
end

%% plot the wheel separating the left and right trials

befWheelAvg = cellfun(@(x)(mean((x),1)),befWheel,'UniformOutput',false);
aftWheelAvg = cellfun(@(x)(mean((x),1)),aftWheel,'UniformOutput',false);

befWheelAvgL = fn_cell2mat(befWheelAvg(:,1),1); befWheelAvgR = fn_cell2mat(befWheelAvg(:,2),1);
aftWheelAvgL = fn_cell2mat(aftWheelAvg(:,1),1); aftWheelAvgR = fn_cell2mat(aftWheelAvg(:,2),1);

figure; 
for i = 1:mouseMega.nMouse
    subplot(2,3,i); 
    p = fn_plotComparison({befWheelAvgL(:,i),aftWheelAvgL(:,i)},'paired',true,'compType','errorbarWithDot');
    xlim([0.6 2.4]); xticks(1:3);title(['p=' num2str(p,'%.4f')]); xticklabels({'bias','unbias'});
end
figure; 
for i = 1:mouseMega.nMouse
    subplot(2,3,i); 
    p = fn_plotComparison({befWheelAvgR(:,i),aftWheelAvgR(:,i)},'paired',true,'compType','errorbarWithDot');
    xlim([0.6 2.4]); xticks(1:3);title(['p=' num2str(p,'%.4f')]); xticklabels({'bias','unbias'});
end

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
