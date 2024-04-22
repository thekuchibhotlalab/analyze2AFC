%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT_bin40(mice);
%for i = 1:length(mice); allMouse{i} = allMouse{i}.getPsyTrack(); end

mouseMega = wheel2AFCmega(allMouse);

%%
accThre = [0.7 nan]; windowSize = 400;tempBehav = cell(1,3);
for i = 1:length(allMouse)
    [~, learningOnset] = allMouse{i}.fn_findLearningOnset(accThre,windowSize);
    selTrial = learningOnset-300:learningOnset+300;
    nTrials = size(allMouse{i}.behav,1); 
    selTrial(selTrial<=0) = []; selTrial(selTrial>nTrials) = [];
    dayStart = allMouse{i}.behav.day(selTrial(1)); dayEnd = allMouse{i}.behav.day(selTrial(end));
    trialLim(1) = find(allMouse{i}.behav.day==dayStart,1); trialLim(2) = find(allMouse{i}.behav.day==dayEnd,1,'last');

    %figure; subplot(2,1,1);plot(allMouse{i}.behav.modelPred_removeAAcc); hold on; plot(allMouse{i}.behav.modelPredAcc); plot(allMouse{i}.behav.acc)
    %plot([trialLim(1) trialLim(1)],[0 1]); plot([trialLim(2) trialLim(2)],[0 1])
    %subplot(2,1,2); plot(abs(allMouse{i}.behav.bias)); hold on; plot(abs(allMouse{i}.behav.modelBias));
    %plot([trialLim(1) trialLim(1)],[0 1]); plot([trialLim(2) trialLim(2)],[0 1])
    
    binSize = 40;
    sampleInterval = 20;
    tempBin = trialLim(1)+sampleInterval:sampleInterval:(trialLim(2)-sampleInterval);
    tempBehav{1}{i} = nan(1,length(tempBin)); tempBehav{2}{i} = nan(1,length(tempBin)); 
    tempBehav{3}{i} = nan(1,length(tempBin)); tempBehav{4}{i} = nan(1,length(tempBin));
    tempStim = [];
    for j = 1:length(tempBin)        
        tempIdx = tempBin(j)-sampleInterval:tempBin(j)+sampleInterval-1;

        [tempBehav{1}{i}(j),tempBehav{2}{i}(j)] = fn_getAccBias(allMouse{i}.behav.stimulus(tempIdx),allMouse{i}.behav.responseType(tempIdx)==1,...
            allMouse{i}.behav.action(tempIdx)==0);
        tempStim = nanmean(allMouse{i}.behav.stimulus(tempIdx)==1);
        tempBehav{3}{i}(j) = allMouse{i}.behav.modelPred_removeACorr(tempBin(j));
        tempBehav{4}{i}(j) = allMouse{i}.behav.modelPred_removeHAcc(tempBin(j));
    end
end 
%%
tempBehavMat = cellfun(@(x)(fn_cell2matFillNan(x)),tempBehav,'UniformOutput',false);
tempBehavAvg = cellfun(@(x)(nanmean(x,2)),tempBehavMat,'UniformOutput',false);
%tempDiff1 = cellfun(@(x,y)(x-y),tempBehav{3}, tempBehav{2},'UniformOutput',false); 

figure; fn_plotComparison({tempBehavAvg{3}-tempBehavAvg{2},tempBehavAvg{4}-tempBehavAvg{2}},'paired',true,'lineType','none','dotType','random','barType','bar',...
    'errorbarArgIn', {'Color',[0.2 0.2 0.2],'LineWidth',1,'LineStyle','none'},'barplotArgIn', {0.4,'EdgeColor','none','FaceColor',matlabColors(2,0.6),'LineWidth',0.5},...
    'scatterArgIn', {15, [0.6 0.6 0.6]}); 
xlim([0.2 2.8]); ylim([-0.05 0.15])
[p1] = signrank(tempBehavAvg{3}-tempBehavAvg{2});  [p2] = signrank(tempBehavAvg{4}-tempBehavAvg{2}); title(['p1=' num2str(p1) ' p2=' num2str(p2)])

figure; hold on;
tempA = abs(tempBehavMat{1}(:)); tempB = tempBehavMat{3}(:) - tempBehavMat{2}(:);
[lm, sc] = fn_plotScatterCorr({tempA,tempA},{tempBehavMat{4}(:) - tempBehavMat{2}(:),tempB},'scatterColors', {[0.5 0.5 0.5],matlabColors(2)}, ...
    'evalCorr', true,'scatterArgIn',{'filled','MarkerFaceAlpha',0.6},'corrPlotArgIn',{'LineWidth',0.8});
xlim([0 1]);

%%
function acc = getBias(tempStim,tempCorrect)
    acc_s1 = nanmean(tempCorrect (tempStim==1) );
    acc_s2 = nanmean(tempCorrect (tempStim==2) );
    bias = (acc_s1 - acc_s2); % / (acc_s1 + acc_s2);
    acc = (acc_s1 + acc_s2)/2;
end