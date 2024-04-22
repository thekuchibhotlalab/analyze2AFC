%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT_bin30(mice);
mouseMega = wheel2AFCmega(allMouse);
% align data by minimizing performance variance
probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
outCell = mouseMega.objFun('binProbeByTrialFromLearningOnset',{[reinfThre probeThre],probeTrialBin});
reinfAlignPoint = cell2mat(outCell{3});
%% FIND THE LOCATION OF THE LONGEST BIAS BLOCK

blockL = getProp(mouseMega,'biasBlock','field','blockL');
blockR = getProp(mouseMega,'biasBlock','field','blockR');

allBlockLen = cellfun(@(x,y)([x.len y.len]),blockL,blockR,'UniformOutput',false);
allBlockStart = cellfun(@(x,y)([x.start y.start]),blockL,blockR,'UniformOutput',false);
allBlockEnd = cellfun(@(x,y)([x.end y.end]),blockL,blockR,'UniformOutput',false);
allBlockAcc = cellfun(@(x,y)([x.acc y.acc]),blockL,blockR,'UniformOutput',false);

longestBlock = [];
for i = 1:length(allBlockLen)
    [tempLen,tempIdx] = max(allBlockLen{i}); 
    longestBlock(i,:) = [allBlockStart{i}(tempIdx)-reinfAlignPoint(i) allBlockEnd{i}(tempIdx)-reinfAlignPoint(i) tempLen allBlockAcc{i}(tempIdx)];
end
longestBlock = array2table(longestBlock,'VariableNames',{'start','end','len','acc'});


%% visualize accuracy and bias curve for biased blocks
mouse = 8; 
figure; 
for i = 1:length(allBlockLen{mouse})
    [temp1,temp2] = fn_sqrtInt(length(allBlockLen{mouse}));
    subplot_tight(temp1,temp2,i,[0.01 0.01]); hold on; 
    plot([0.2 1.0],[0 0],'Color',[0.8 0.8 0.8]);plot([0 0],[-1 1],'Color',[0.8 0.8 0.8])
    tempTrial = [allBlockStart{mouse}(i) allBlockEnd{mouse}(i)];
    temp = getSubsampleAccBias(allMouse,mouse,tempTrial); plot(temp{:},'-o','MarkerSize',3);
    xlim([0.2 1.0]); ylim([-1 1])
end
%% visualize accuracy and bias curve for trials in training
mouse = 13; binSize = 0:200:size(allMouse{mouse}.behav,1);
figure; tempMean = []; tempSEM = [];
[temp1,temp2] = fn_sqrtInt(length(binSize));
for i = 1:length(binSize)-1
    subplot_tight(temp1,temp2,i,[0.01 0.01]); hold on; 
    plot([0.2 1.0],[0 0],'Color',[0.8 0.8 0.8]);plot([0 0],[-1 1],'Color',[0.8 0.8 0.8])
    tempTrial = [binSize(i)+1 binSize(i+1)];
    [temp,tempMean(i), tempSEM(i)] = getSubsampleAccBias(allMouse,mouse,tempTrial); plot(temp{:},'-o','MarkerSize',3);
    xlim([0.2 1.0]); ylim([-1 1])
    legend({'','',['acc=' num2str(nanmean(allMouse{mouse}.behav.acc(tempTrial(1):tempTrial(2))))]})
    xticks([]); yticks([]);
end
%subplot_tight(temp1,temp2,i+1,[0.01 0.01]);
%plot(tempMean)

%% visualize random choice trajectory 
randStim = rand(1,2000); binSize = 0:200:2000; 
randStim(randStim<0.5) = -1; randStim(randStim>=0.5) = 1;
logW = linspace(0,0,2000);
choiceProb = 1./(1+exp(-logW.*randStim));
randStim1 = randStim == -1; randStim2 = randStim == 1;

randChoice = rand(1,2000);c2 = (randChoice<=choiceProb); c1 = (randChoice>choiceProb);
randChoice(c1) = -1; randChoice(c2) = 1;
acc1 = (randChoice==randStim).*randStim1; acc1(randStim2) = nan;
acc2 =  (randChoice==randStim).*randStim2; acc2(randStim1) = nan;
acc1 = smoothdata(acc1,'gaussian',30); acc2 = smoothdata(acc2,'gaussian',30);
accracy = acc1/2+acc2/2; bias = smoothdata(acc2-acc1,'gaussian',30);
figure; tempMean = []; tempSEM = [];
for i = 1:length(binSize)-1
    [temp1,temp2] = fn_sqrtInt(length(binSize));
    subplot_tight(temp1,temp2,i,[0.01 0.01]); hold on; 
    plot([0.2 1.0],[0 0],'Color',[0.8 0.8 0.8]);plot([0 0],[-1 1],'Color',[0.8 0.8 0.8])
    tempTrial = [binSize(i)+1 binSize(i+1)];
    plot(accracy(tempTrial(1):tempTrial(2)),bias(tempTrial(1):tempTrial(2)),'-o','MarkerSize',3);

    tempXDiff = diff(accracy(tempTrial(1):tempTrial(2))); tempYDiff = diff(bias(tempTrial(1):tempTrial(2)));
    spd = sqrt(tempXDiff.^2 + tempYDiff.^2);
    tempMean(i) = nanmean(spd); tempSEM(i) = nanstd(spd) ./ sqrt(length(spd)); 

    xlim([0.2 1.0]); ylim([-1 1])
    legend({'','',['acc=' num2str(nanmean(accracy(tempTrial(1):tempTrial(2))))]})
    xticks([]); yticks([]);
end
subplot_tight(temp1,temp2,i+1,[0.01 0.01]);
plot(tempMean)

%% Visualize acc1 and acc2 curve

%{
figure; tempSize = size(allMouse{mouse}.biasBlock.trans,1);
for i = 1:tempSize
    [temp1,temp2] = fn_sqrtInt(tempSize);
    subplot_tight(temp1,temp2,i,[0.01 0.01]); hold on; 
    plot([0 1.0],[0.5 0.5],'Color',[0.8 0.8 0.8]);plot([0.5 0.5],[0 1],'Color',[0.8 0.8 0.8])
    tempTrial = [allMouse{mouse}.biasBlock.trans(i,1) allMouse{mouse}.biasBlock.trans(i,2)];
    temp = getSubsampleAccLR(allMouse,mouse,tempTrial); plot(temp{:},'-o','MarkerSize',3);
    xlim([0 1.0]); ylim([0 1])
end
%}

%% ALL FUNCTIONS

function [temp,spdmean, spdSEM] = getSubsampleAccBias(allMouse,mouse,trialLim)
temp = {};
trialLim(1) = max([trialLim(1) 1]);
trialLim(2) = min([trialLim(2) length(allMouse{mouse}.behav.bias)]);
sampleBin = 4; sampleAxis = trialLim(1):sampleBin:trialLim(2);
tempX = smoothdata(allMouse{mouse}.behav.acc,'gaussian',20);
tempY = smoothdata(allMouse{mouse}.behav.bias,'gaussian',20);
temp{1} = tempX(sampleAxis); temp{2} = tempY(sampleAxis);
tempXDiff = diff(tempX(sampleAxis)); tempYDiff = diff(tempY(sampleAxis));
spd = sqrt(tempXDiff.^2 + tempYDiff.^2);
spdmean = nanmean(spd); spdSEM = nanstd(spd) ./ sqrt(length(spd)); 
end
function temp = getSubsampleAccLR(allMouse,mouse,trialLim)
temp = {};
trialLim(1) = max([trialLim(1) 1]);
trialLim(2) = min([trialLim(2) length(allMouse{mouse}.behav.bias)]);
sampleBin = 4; sampleAxis = trialLim(1):sampleBin:trialLim(2);
tempX = smoothdata(allMouse{mouse}.behav.acc1,'gaussian',20);
tempY = smoothdata(allMouse{mouse}.behav.acc2,'gaussian',20);
temp{1} = tempX(sampleAxis); temp{2} = tempY(sampleAxis);
end





