%% LOAD DATA
clear; 
mice ={'zz107','zz109','zz111','zz112','zz113','zz115'};
opsParam.biasBlockType = 'threshold';
allMouse = fn_getObjPT_bin40(mice,opsParam); mouseMega = wheel2AFCmega(allMouse);

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
for i = 1:mouseMega.nMouse
    accThre = [0.7 nan]; windowSize= 400; 
    [~, learningOnset] = allMouse{i}.fn_findLearningOnset(accThre,windowSize);
    selTrial = learningOnset-300:learningOnset+300;
    nTrials = size(allMouse{i}.behav,1); 
    selTrial(selTrial<=0) = []; selTrial(selTrial>nTrials) = [];
    dayStart = allMouse{i}.behav.day(selTrial(1)); dayEnd = allMouse{i}.behav.day(selTrial(end));
    
    tempSelIdx = false(nTrials,1);
    tempSelIdx(~badFlag{i} & allMouse{i}.behav.day>=dayStart & allMouse{i}.behav.day<=dayEnd & (allMouse{i}.behav.reinfBef | allMouse{i}.behav.reinfAft | allMouse{i}.behav.probe)) = true;

    tempSpeed = diff(allMouse{i}.behavVar.wheelPos_aligned,1,2);
    preMoveAbs = nanmean(abs(tempSpeed(tempSelIdx,1:3000)),2);
    preMove = nanmean((tempSpeed(tempSelIdx,1:3000)),2);
    allTbl_ctxt{i} = allMouse{i}.behav(tempSelIdx,:);
    allTbl_ctxt{i}.label = labelCell{i}(tempSelIdx);
    for j = 1:length(wheelAttrName); allTbl_ctxt{i}.(wheelAttrName{j}) = attributesCell{i}(tempSelIdx,j); end

    allTbl_ctxt{i}.preMoveAbs = preMoveAbs;
    allTbl_ctxt{i}.preMove = preMove .* (allTbl_ctxt{i}.action*2-3);
end
reinfBefRT = zeros(mouseMega.nMouse,2); reinfAftRT = zeros(mouseMega.nMouse,2); probeRT = zeros(mouseMega.nMouse,2);
reinfBefType = zeros(mouseMega.nMouse,1); reinfAftType = zeros(mouseMega.nMouse,1); probeType = zeros(mouseMega.nMouse,1);
reinfBefMov = zeros(mouseMega.nMouse,1); reinfAftMov = zeros(mouseMega.nMouse,1); probeMov = zeros(mouseMega.nMouse,1);

[rpRT] = getTableItem(allTbl_ctxt,'reactionTime',true);
[rpRT_all] = getTableItem(allTbl_ctxt,'reactionTime',false);
[rpLabel] = getTableItem(allTbl_ctxt,'label',false);
[rpPreMove] = getTableItem(allTbl_ctxt,'preMove',false);
[rpITIMove] = getTableItem(allTbl_ctxt,'itiMove',false);
[rpITIChoice] = getTableItem(allTbl_ctxt,'itiChoiceDir',false);
[rpOverShoot] = getTableItem(allTbl_ctxt,'overshoot',false);

rpWheel = {}; 
for i = 1:length(wheelAttrName); rpWheel{i} = getTableItem(allTbl_ctxt,wheelAttrName{i},false); end

%{
figure; subplot(2,4,1)
p = fn_plotComparison({rpRT.bef(:,1),rpRT.probe(:,1),rpRT.aft(:,1)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 3.4]); xticks(1:3); ylim([0.3 0.7]); 
subplot(2,4,2)
p = fn_plotComparison({rpRT_all.bef(:,1),rpRT_all.probe(:,1),rpRT_all.aft(:,1)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 3.4]); xticks(1:3)

for i = 1:length(rpWheel)
    subplot(2,4,i+2);
    p = fn_plotComparison({rpWheel{i}.bef,rpWheel{i}.probe,rpWheel{i}.aft},'paired',true,'compType','errorbarWithDot');
    xlim([0.6 3.4]); xticks(1:3)
end
%}

figure; subplot(2,4,1)
p = fn_plotComparison({rpRT.bef(:,1),rpRT.probe(:,1)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks(1:2); ylim([0.3 0.7]); yticks(0.3:0.2:0.7)
subplot(2,4,2)
p = fn_plotComparison({rpRT_all.bef(:,1),rpRT_all.probe(:,1)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks(1:2); ylim([0.3 1.1]); yticks(0.3:0.4:1.1)

ylimCell = {[0 800],[3 6],[0.8 1.2],[0 500],[3 10],[1 1.5]};
for i = 1:length(rpWheel)
    subplot(2,4,i+2);
    p = fn_plotComparison({rpWheel{i}.bef,rpWheel{i}.probe},'paired',true,'compType','errorbarWithDot');
    xlim([0.6 2.4]); xticks(1:2); ylim(ylimCell{i}); 
end
%%
figure; subplot(1,4,1);
p = fn_plotComparison({rpLabel.bef(:,1),rpLabel.probe(:,1)},'paired',true,'compType','errorbarWithDot');
subplot(1,4,2);
p = fn_plotComparison({rpITIChoice.bef(:,1),rpITIChoice.probe(:,1)},'paired',true,'compType','errorbarWithDot');
subplot(1,4,3);
p = fn_plotComparison({rpOverShoot.bef(:,1),rpOverShoot.probe(:,1)},'paired',true,'compType','errorbarWithDot');

%%
figure;
p = fn_plotComparison({rpPreMove.bef,rpPreMove.probe,rpPreMove.aft},'paired',true,'compType','errorbarWithDot');
xlim([0.6 3.4]); xticks(1:3)
%% by session
reinfRT = cell(mouseMega.nMouse,1); probeRT = cell(mouseMega.nMouse,1);
reinfType = cell(mouseMega.nMouse,1); probeType = cell(mouseMega.nMouse,1);
reinfPerf = cell(mouseMega.nMouse,1); probePerf = cell(mouseMega.nMouse,1);
for i = 1:mouseMega.nMouse
    days = unique(allTbl_ctxt{i}.day);
    tempProbeRTDay = zeros(length(days),2);tempReinfRTDay = zeros(length(days),2);
    tempProbePerfDay = zeros(length(days),2);tempReinfPerfDay = zeros(length(days),2);
    for day = 1:length(days)
        
        tempProbeRT = nan(2,2); tempReinfRT = nan(2,2);
        tempBehav = allTbl_ctxt{i}(allTbl_ctxt{i}.day==days(day),:);
        for dir = 1:2 % loop through action type          
            for l = 1:2 % loop through label
                tempIdx =  tempBehav.action == dir & tempBehav.label == (l-1);
                tempProbeRT(l,dir) =  nanmean(tempBehav.reactionTime(tempIdx & tempBehav.goodProbe));
                tempReinfRT(l,dir) =  nanmean(tempBehav.reactionTime(tempIdx & tempBehav.reinfBef));
            end 
        end
        tempReinfRTDay(day,:)= nanmean(tempReinfRT,2);
        tempProbeRTDay(day,:) = nanmean(tempProbeRT,2);

        tempBehavR = tempBehav(logical(tempBehav.reinfBef),:);
        [tempReinfPerfDay(day,2),tempReinfPerfDay(day,1)] = fn_getAccBias(tempBehavR.stimulus,tempBehavR.responseType==1);
        tempBehavP = tempBehav(logical(tempBehav.goodProbe),:);
        [tempProbePerfDay(day,2),tempProbePerfDay(day,1)] = fn_getAccBias(tempBehavP.stimulus,tempBehavP.responseType==1);
        reinfType{i}(day,1) = mean(tempBehavR.label);
        probeType{i}(day,1) = mean(tempBehavP.label);
    end
    reinfRT{i} = tempReinfRTDay; probeRT{i} = tempProbeRTDay;
    reinfPerf{i} = tempReinfPerfDay; probePerf{i} = tempProbePerfDay;
end
reinfRT = fn_cell2mat(reinfRT,1); probeRT = fn_cell2mat(probeRT,1);
reinfPerf = fn_cell2mat(reinfPerf,1); probePerf = fn_cell2mat(probePerf,1);
reinfType = fn_cell2mat(reinfType,1); probeType = fn_cell2mat(probeType,1);

a = reinfRT-probeRT;
b = probePerf-reinfPerf;
c = reinfType-probeType;
figure; scatter(a(:,1),b(:,1));
figure; fn_plotComparison({reinfRT(:,2),probeRT(:,2)},'paired',true);
figure; fn_plotComparison({reinfType,probeType},'paired',true);
%% FUNCTIONS
function [item] = getTableItem(allTbl_ctxt,itemName,separateLabel)
if nargin == 2; separateLabel = false; end
if separateLabel
    item.bef = zeros(length(allTbl_ctxt),2); item.aft = zeros(length(allTbl_ctxt),2); item.probe = zeros(length(allTbl_ctxt),2);
else
    item.bef = zeros(length(allTbl_ctxt),1); item.aft = zeros(length(allTbl_ctxt),1); item.probe = zeros(length(allTbl_ctxt),1);
end
for i = 1:length(allTbl_ctxt)
    if separateLabel
        tempBef = nan(2,2); tempAft = nan(2,2); tempProbe = nan(2,2);
    else
        tempBef = nan(2,1); tempAft = nan(2,1);tempProbe = nan(2,1);
    end

    for dir = 1:2 % loop through action type      
        if separateLabel
            for l = 1:2 % loop through label
                tempIdx =  allTbl_ctxt{i}.action == dir & allTbl_ctxt{i}.label == (l-1);
                tempProbe(dir,l) =  nanmean(allTbl_ctxt{i}.(itemName)(tempIdx & allTbl_ctxt{i}.goodProbe));
                tempBef(dir,l) =  nanmean(allTbl_ctxt{i}.(itemName)(tempIdx & allTbl_ctxt{i}.reinfBef));
                tempAft(dir,l) =  nanmean(allTbl_ctxt{i}.(itemName)(tempIdx & allTbl_ctxt{i}.reinfAft));
            end
        else
            tempBef(dir) = mean(allTbl_ctxt{i}.(itemName)(allTbl_ctxt{i}.action == dir & allTbl_ctxt{i}.reinfBef));
            tempProbe(dir) = mean(allTbl_ctxt{i}.(itemName)(allTbl_ctxt{i}.action == dir & allTbl_ctxt{i}.goodProbe));
            tempAft(dir) = mean(allTbl_ctxt{i}.(itemName)(allTbl_ctxt{i}.action == dir & allTbl_ctxt{i}.reinfAft));
        end  
    end
    item.bef(i,:)= nanmean(tempBef,1); item.aft(i,:)= nanmean(tempAft,1); item.probe(i,:) = nanmean(tempProbe,1); 
end

end
