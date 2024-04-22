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

[rpRT_all] = getTableItem(allTbl_ctxt,'reactionTime',false);
[rpLabel] = getTableItem(allTbl_ctxt,'label',false);
[rpPreMoveAbs] = getTableItem(allTbl_ctxt,'preMoveAbs',false);
[rpITIMoveAbs] = getTableItem(allTbl_ctxt,'itiMoveAbs',false);
[rpITIChoice] = getTableItem(allTbl_ctxt,'itiChoiceDir',false);
[rpOverShoot] = getTableItem(allTbl_ctxt,'overshoot',false);

rpWheel = {}; 
for i = 1:length(wheelAttrName); rpWheel{i} = getTableItem(allTbl_ctxt,wheelAttrName{i},false); end



%%

figure; subplot(1,4,1); hold on;
fn_plotComparison( (rpRT_all.probe./rpRT_all.bef-1)','barType','bar','barplotArgIn', ...
    {0.4,'EdgeColor','none','FaceColor',matlabColors(2,0.6),'LineWidth',0.5},'lineType',{},'errorbarArgIn',...
    {'Color',[0.2 0.2 0.2],'LineWidth',0.5},'scatterArgIn',{15, [0.6 0.6 0.6],'filled'},'dotType','random'); 
xlim([0 2]); ylim([-0.4 0.4]);title('reaction time')

subplot(1,4,2)
fn_plotComparison((rpLabel.probe-rpLabel.bef)','barType','bar','barplotArgIn', ...
    {0.4,'EdgeColor','none','FaceColor',matlabColors(2,0.6),'LineWidth',0.5},'lineType',{},'errorbarArgIn',...
    {'Color',[0.2 0.2 0.2],'LineWidth',0.5},'scatterArgIn',{15, [0.6 0.6 0.6],'filled'},'dotType','random');
xlim([0 2]); ylim([-0.2 0.2]); title('label')


subplot(1,4,3)
fn_plotComparison((rpWheel{5}.probe./rpWheel{5}.bef-1)','barType','bar','barplotArgIn', ...
    {0.4,'EdgeColor','none','FaceColor',matlabColors(2,0.6),'LineWidth',0.5},'lineType',{},'errorbarArgIn',...
    {'Color',[0.2 0.2 0.2],'LineWidth',0.5},'scatterArgIn',{15, [0.6 0.6 0.6],'filled'},'dotType','random');
xlim([0 2]); ylim([-0.5 0.5]); title('total speed')

subplot(1,4,4)
fn_plotComparison((rpWheel{2}.probe./rpWheel{2}.bef-1)','barType','bar','barplotArgIn', ...
    {0.4,'EdgeColor','none','FaceColor',matlabColors(2,0.6),'LineWidth',0.5},'lineType',{},'errorbarArgIn',...
    {'Color',[0.2 0.2 0.2],'LineWidth',0.5},'scatterArgIn',{15, [0.6 0.6 0.6],'filled'},'dotType','random');
xlim([0 2]); ylim([-0.5 0.5]); title('onset speed')

%% 
figure; subplot(2,1,1)
fn_plotComparison((rpRT_all.probe./rpRT_all.bef-1)','paired',false,'lineType','none','dotType','random','barType','bar','compType','errorbarWithDot',...
    'errorbarArgIn', {'Color',[0.2 0.2 0.2],'LineWidth',1,'LineStyle','none'},'barplotArgIn', {0.4,'EdgeColor','none','FaceColor',matlabColors(2,0.6),'LineWidth',0.5},...
    'scatterArgIn', {15, [0.6 0.6 0.6]}); 
xlim([0.2 1.8]); ylim([-0.3 0.3]); ylabel('reaction time'); [h,p] = ttest(rpRT_all.probe./rpRT_all.bef-1);

subplot(2,1,2)
fn_plotComparison((rpWheel{5}.probe./rpWheel{5}.bef-1)','paired',false,'lineType','none','dotType','random','barType','bar','compType','errorbarWithDot',...
    'errorbarArgIn', {'Color',[0.2 0.2 0.2],'LineWidth',1,'LineStyle','none'},'barplotArgIn', {0.4,'EdgeColor','none','FaceColor',matlabColors(2,0.6),'LineWidth',0.5},...
    'scatterArgIn', {15, [0.6 0.6 0.6]}); ylim([-0.1 0.4]); ylabel('total speed')
xlim([0.2 1.8]); p = signrank(rpWheel{5}.probe./rpWheel{5}.bef-1);
%% function

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
