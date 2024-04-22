clear;
global dataPath sep;
dataPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\trialData';
sep = '\';
%selectProtocol = {'FM_One_Oct','FM_One_Oct_Prob','FM_Half_Oct','FM_Half_Oct_Prob'};
selectProtocol = {'puretone'};

%% PLOT ACCURACY LEARNING CURVE ACROSS MICE
%mice = {'zz062','zz063'};
mice = {'zz054','zz062','zz063','zz066','zz067','zz068','zz069'};
%mice = {'zz071','zz073','zz075','zz076','zz077','zz048','zz050'};
%mice = {'zz070','zz071','zz072','zz073'};

[accuracy,accuracyTop,accuracyMean,acc_L,accMean_L,acc_R,accMean_R,...
    bias,bias_absMean,missBias,probeData,reinfDataBef, reinfDataAft,...
    trialNum,probeTrialNum,probeDayNum,dayLen] = ...
    cellfun(@(x)(getDataByAnimal(x,selectProtocol, true)),mice,'UniformOutput',false);

bias_abs = cellfun(@abs,bias,'UniformOutput',false);
%{
for i = 1:length(mice)
    mouse = mice{i};
    load([dataPath sep mouse '.mat']);
    
    responseType = trialData.responseType; 

    
    [actionRate, ~] = fn_cell2mat(cellfun(@(x)(smoothdata(double(x~=0),'movmean',...
        learningCurveBin)),responseType,'UniformOutput',false),1);
    [accuracy, cellSize] = fn_cell2mat(cellfun(@(x)(smoothdata(double(x==1),'movmean',...
        learningCurveBin)),responseType,'UniformOutput',false),1);
    accuracyCell{i} = accuracy./actionRate;
    
    [bias, acc_L, acc_R] = cellfun(@getBias,trialData.stimulus,trialData.responseType,'UniformOutput',false);
    bias = fn_cell2mat(bias,1); acc_L = fn_cell2mat(acc_L,1); acc_R = fn_cell2mat(acc_R,1);
    biasCell{i} = abs(bias); accCell_L{i} = acc_L; accCell_R{i} = acc_R;
    
    missBias = cellfun(@getMissBias,trialData.stimulus,trialData.responseType,'UniformOutput',false);
    missBias = fn_cell2mat(missBias,1);
    %--------PROBE ORGANIZE BY DAY---------
    [probeData, trialNum,probeTrialNum,trialNumNoMiss,probeTrialNumNoMiss] = cellfun(@getProbe,trialData.stimulus,trialData.responseType,...
        trialData.context,'UniformOutput',false);
    probeData = fn_cell2mat(probeData,1);
    trialNum = fn_cell2mat(trialNum,1);
    probeTrialNum = fn_cell2mat(probeTrialNum,1);
    
    cumsumTrials = [0;cumsum(trialNum)];
    probeTrialNum = probeTrialNum + cumsumTrials(1:end-1);
    probeTrialNum(isnan(probeTrialNum)) = [];
    %-------------------------------
    
    %stimulus = fn_cell2mat(trialData.stimulus);
    %responseType = fn_cell2mat(trialData.responseType);
    %context = fn_cell2mat(trialData.context);
    %probeFlag = context==3; probeStart = find(diff(probeFlag)==1)+1; probeEnd = find(diff(probeFlag)==-1);
    %nProbe = probeEnd-probeStart+1; extraTrialFlag = double(mod(nProbe,10)==1); 
    %probeStart = probeStart + extraTrialFlag;
    %probePerformance = getProbePerformance(stimulus,responseType,probeStart,probeEnd);
    
    probeCell{i} = probeData; %probePerformance;
    probeTrialNumCell{i} = probeTrialNum;%(probeStart+probeEnd)/2;
    
    dayLen = cellSize(1,:); dayLen(dayLen==0) = [];
    
    plotAnimalByDay([accuracy actionRate bias missBias],probeData,probeTrialNum,dayLen,...
        {'Accuracy','ActionRate','ActionBias','MissBias'},mouse)
end
%}
%% PLOT ALL ANIMALS PERORMANCE
plotAllAnimalsReinfProbe(accuracy,probeData,probeTrialNum);
%% PLOT ALL ANIMALS BIAS
plotAllAnimalsReinfProbeBias(bias_abs,probeData,probeTrialNum);

%% PLOT ALL ANIMALS PERFORMANCE BY DAY
nDay = 8;
accuracyMeanPlot = cellfun(@fn_removeNan,accuracyMean,'UniformOutput',false);
accuracyMeanPlot = fn_cell2mat(cellfun(@(x)(x(1:nDay)),accuracyMeanPlot,'UniformOutput',false),2);

accuracyTopPlot = cellfun(@fn_removeNan,accuracyTop,'UniformOutput',false);
accuracyTopPlot = fn_cell2mat(cellfun(@(x)(x(1:nDay)),accuracyTopPlot,'UniformOutput',false),2);

probeDayIdx = cellfun(@(x)(x<=nDay),probeDayNum,'UniformOutput',false);
probeDayAlign = nan(nDay,length(probeData),size(probeData{1},2));
reinfBefDayAlign = nan(nDay,length(probeData),size(probeData{1},2));
reinfAftDayAlign = nan(nDay,length(probeData),size(probeData{1},2));


for i = 1:length(probeData)
    probeDayAlign(probeDayNum{i}(probeDayIdx{i}),i,:) = probeData{i}(probeDayIdx{i},:);
    reinfBefDayAlign(probeDayNum{i}(probeDayIdx{i}),i,:) = reinfDataBef{i}(probeDayIdx{i},:);
    reinfAftDayAlign(probeDayNum{i}(probeDayIdx{i}),i,:) = reinfDataAft{i}(probeDayIdx{i},:);
end
probeDayPlot = mean(probeDayAlign(:,:,1:2),3);
reinfBefDayPlot = mean(reinfBefDayAlign(:,:,1:2),3);
reinfAftDayPlot = mean(reinfAftDayAlign(:,:,1:2),3);


% PERFORMANCE OVER DAYS
figure; 
subplot(2,2,1);hold on;
f1 = fn_plotMeanErrorbar(1:nDay,accuracyMeanPlot',matlabColors(1),...
    {'Color',matlabColors(1),'LineWidth',2},{'faceAlpha',0.2,'LineStyle','none'});
f2 = fn_plotMeanErrorbar(1:nDay,accuracyTopPlot',matlabColors(6),...
    {'Color',matlabColors(6),'LineWidth',2},{'faceAlpha',0.2,'LineStyle','none'});
f3 = fn_plotMeanErrorbar(1:nDay,probeDayPlot',matlabColors(2),...
    {'--o','Color',matlabColors(2),'LineWidth',2,},{'faceAlpha',0.2,'LineStyle','none'});
xlabel('Days'); ylabel('Performance'); ylim([0.4 1]); xlim([1 nDay])
legend([f1 f2 f3],'mean reinf','top 10% reinf','probe','Location','Best');

% BAR PLOT
nDay = 7;
subplot(2,2,2); hold on;
tempReinfFlat = reshape(accuracyMeanPlot(1:nDay,:),1,[]);
tempProbeFlat = reshape(probeDayPlot(1:nDay,:),1,[]);
nanFlag = isnan(tempProbeFlat); tempReinfFlat(nanFlag) = []; tempProbeFlat(nanFlag) = [];
[h,p] = ttest(tempReinfFlat,tempProbeFlat,'tail','left');
bar([nanmean(tempReinfFlat) nanmean(tempProbeFlat)] ,'EdgeColor','None','FaceColor',matlabColors(2,0.9)); hold on;
for i = 1:length(tempReinfFlat)
    f = plot([1 2],[tempReinfFlat(i) tempProbeFlat(i)],'Color',[0.6 0.6 0.6],'Marker','.','MarkerSize',15,...
        'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',0.5);
end
legend(f,['p = ' num2str(p,'%.2e')],'Location','Best')
ylabel('Accuracy'); xticks([1 2]); xticklabels({'Reinf','Probe'}); xlim([0 3]);

% SCATTER PLOT
tempColor = reshape(matlabColors(1:size(accuracyMeanPlot,2)),[1 size(accuracyMeanPlot,2) 3]);
tempColor = repmat(tempColor,[nDay 1 1]);
tempColor = reshape(tempColor, [nDay * size(accuracyMeanPlot,2) 3]);
tempColor(nanFlag,:) = [];
subplot(2,2,3);hold on;
limm = [0.4 1];
plot(limm,limm,'Color',[0.8 0.8 0.8],'LineWidth',2);
scatter(tempReinfFlat,tempProbeFlat,30,tempColor,'filled')
xlim(limm); ylim(limm); xlabel('Mean Reinf Accuracy'); ylabel('Mean Probe Accuracy')

tempReinfFlatTop = reshape(accuracyTopPlot(1:nDay,:),1,[]);
tempReinfFlatTop(nanFlag) = [];
subplot(2,2,4);hold on; 
limm = [0.55 1];
plot(limm,limm,'Color',[0.8 0.8 0.8],'LineWidth',2);
scatter(tempReinfFlatTop,tempProbeFlat,30,tempColor,'filled')
xlim(limm); ylim(limm); xlabel('Top 10% Reinf Accuracy'); ylabel('Mean Probe Accuracy')

% BEFORE AND AFTER PLOT
figure; subplot(1,3,1); hold on;
tempBefFlat = reshape(reinfBefDayPlot(1:nDay,:),1,[]);
tempAftFlat = reshape(reinfAftDayPlot(1:nDay,:),1,[]);
tempProbeFlat = reshape(probeDayPlot(1:nDay,:),1,[]);
nanFlag = isnan(tempProbeFlat); 
tempBefFlat(nanFlag) = []; tempAftFlat(nanFlag) = []; tempProbeFlat(nanFlag) = [];
[hBef,pBef] = ttest(tempBefFlat,tempProbeFlat,'tail','left');
[hAft,pAft] = ttest(tempAftFlat,tempProbeFlat,'tail','left');
bar([nanmean(tempBefFlat) nanmean(tempProbeFlat) nanmean(tempAftFlat)] ,'EdgeColor','None','FaceColor',matlabColors(2,0.9)); hold on;
for i = 1:length(tempBefFlat)
    f = plot([1 2 3],[tempBefFlat(i) tempProbeFlat(i) tempAftFlat(i)],'Color',[0.6 0.6 0.6],'Marker','.','MarkerSize',15,...
        'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',0.5);
end
legend(f,['pBef = ' num2str(pBef,'%.2e') newline 'pAft = ' num2str(pAft,'%.2e')],'Location','Best')
ylabel('Accuracy'); xticks([1 2 3]); xticklabels({'Bef','Probe','Aft'}); xlim([0 4]);

subplot(1,3,2); hold on;
limm = [0.4 1];
plot(limm,limm,'Color',[0.8 0.8 0.8],'LineWidth',2);
scatter(tempBefFlat,tempProbeFlat,30,tempColor,'filled');
xlim(limm); ylim(limm); xlabel('Reinf Accuracy BEF'); ylabel('Mean Probe Accuracy')

subplot(1,3,3); hold on;
limm = [0.4 1];
plot(limm,limm,'Color',[0.8 0.8 0.8],'LineWidth',2);
scatter(tempAftFlat,tempProbeFlat,30,tempColor,'filled');
xlim(limm); ylim(limm); xlabel('Reinf Accuracy AFT'); ylabel('Mean Probe Accuracy')
%% PLOT ALL ANIMALS BIAS BY DAY
nDay = 7;
bias_absPlot = fn_cell2mat(cellfun(@(x)(x(1:nDay)),bias_absMean,'UniformOutput',false),2);
probeBiasPlot = abs(probeDayAlign(1:nDay,:,5));
reinfBefBiasPlot = abs(reinfBefDayAlign(1:nDay,:,5));
reinfAftBiasPlot = abs(reinfAftDayAlign(1:nDay,:,5));

tempPerf = accuracyMeanPlot(1:nDay,:); tempSize = tempPerf.^2 *80;
tempColor = reshape(matlabColors(1:size(bias_absPlot,2)),[1 size(bias_absPlot,2) 3]);
tempColor = repmat(tempColor,[size(bias_absPlot,1) 1 1]);
tempColor = reshape(tempColor, [size(bias_absPlot,1) * size(bias_absPlot,2) 3]);

% PLOT BIAS SCATTER PLOT
figure; subplot(2,2,1); hold on
maxLim = round(max([bias_absPlot(:);probeBiasPlot(:)])/0.1)*0.1;
plot([0 maxLim],[0 maxLim],'Color',[0.8 0.8 0.8],'LineWidth',2);
scatter(bias_absPlot(:),probeBiasPlot(:),round(tempSize(:)),tempColor,'filled');
xlim([0 maxLim]); ylim([0 maxLim])
xlabel('Mean Bias Reinf'); ylabel('Mean Bias Probe'); title('Mean Bias (Abs)')

accMeanPlot_L = fn_cell2mat(cellfun(@(x)(x(1:nDay)),accMean_L,'UniformOutput',false),2);
accMeanPlot_R = fn_cell2mat(cellfun(@(x)(x(1:nDay)),accMean_R,'UniformOutput',false),2);

probePlotL = probeDayAlign(1:nDay,:,1);probePlotL(probeDayAlign(1:nDay,:,3)<0.4) = nan;
probePlotR = probeDayAlign(1:nDay,:,2);probePlotR(probeDayAlign(1:nDay,:,4)<0.4) = nan;

tempPerf = accuracyMeanPlot(1:nDay,:); tempSize = tempPerf.^2 *80;
tempColor = reshape(matlabColors(1:size(accMeanPlot_L,2)),[1 size(accMeanPlot_L,2) 3]);
tempColor = repmat(tempColor,[size(accMeanPlot_L,1) 1 1]);
tempColor = reshape(tempColor, [size(accMeanPlot_L,1) * size(accMeanPlot_L,2) 3]);

% BAR PLOT
subplot(2,2,2); hold on;
tempReinfFlat = bias_absPlot(:); tempProbeFlat = probeBiasPlot(:); nanFlag = isnan(tempProbeFlat);
tempReinfFlat(nanFlag) = []; tempProbeFlat(nanFlag) = [];
[h,p] = ttest(tempReinfFlat,tempProbeFlat,'tail','right');
bar([nanmean(tempReinfFlat) nanmean(tempProbeFlat)] ,'EdgeColor','None','FaceColor',matlabColors(2,0.9)); hold on;
for i = 1:length(tempReinfFlat)
    f = plot([1 2],[tempReinfFlat(i) tempProbeFlat(i)],'Color',[0.6 0.6 0.6],'Marker','.','MarkerSize',15,...
        'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',0.5);
end
legend(f,['p = ' num2str(p,'%.4f')])
ylabel('Bias'); xticks([1 2]); xticklabels({'Reinf','Probe'}); xlim([0 3]);title('Mean Bias (Abs)')

% BIAS SCATTER PLOT
subplot(2,2,3);hold on;
plot([-1 1],[-1 1],'--','Color',[0.8 0.8 0.8],'LineWidth',1);
plot([0 0],[-1 1],'Color',[0.8 0.8 0.8],'LineWidth',2);
plot([-1 1],[0 0],'Color',[0.8 0.8 0.8],'LineWidth',2);
scatter(accMeanPlot_L(:) - accMeanPlot_R(:),probePlotL(:)-probePlotR(:),tempSize(:),tempColor,'filled');
xlabel('Mean Bias Reinf'); ylabel('Bias Probe'); title('Mean Bias')

% ACC L/R SCATTER PLOT
subplot(2,2,4); hold on
maxLim = round(max([accMeanPlot_L(:);accMeanPlot_R(:)])/0.1)*0.1;
plot([0 maxLim],[0 maxLim],'Color',[0.8 0.8 0.8],'LineWidth',2);
f1 = scatter(accMeanPlot_L(:),accMeanPlot_R(:),30,matlabColors(1),'filled');
f2 = scatter(probePlotL(:),probePlotR(:),30,matlabColors(2),'filled');
legend([f1,f2],'reinf','probe','Location','Best')
xlabel('Left Accuracy'); ylabel('Right Accuracy'); title('Mean Accuracy')

% BEF AND AFT PLOT
figure; subplot(1,3,1); hold on;
tempBefFlat = reshape(reinfBefBiasPlot(1:nDay,:),1,[]);
tempAftFlat = reshape(reinfAftBiasPlot(1:nDay,:),1,[]);
tempProbeFlat = reshape(probeBiasPlot(1:nDay,:),1,[]);
nanFlag = isnan(tempProbeFlat); 
tempBefFlat(nanFlag) = []; tempAftFlat(nanFlag) = []; tempProbeFlat(nanFlag) = [];
[hBef,pBef] = ttest(tempBefFlat,tempProbeFlat,'tail','right');
[hAft,pAft] = ttest(tempAftFlat,tempProbeFlat,'tail','right');
bar([nanmean(tempBefFlat) nanmean(tempProbeFlat) nanmean(tempAftFlat)] ,'EdgeColor','None','FaceColor',matlabColors(2,0.9)); hold on;
for i = 1:length(tempBefFlat)
    f = plot([1 2 3],[tempBefFlat(i) tempProbeFlat(i) tempAftFlat(i)],'Color',[0.6 0.6 0.6],'Marker','.','MarkerSize',15,...
        'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',0.5);
end
legend(f,['pBef = ' num2str(pBef,'%.2e') newline 'pAft = ' num2str(pAft,'%.2e')],'Location','Best')
ylabel('Bias(Abs)'); xticks([1 2 3]); xticklabels({'Bef','Probe','Aft'}); xlim([0 4]);

subplot(1,3,2); hold on
maxLim = round(max([tempBefFlat(:);tempBefFlat(:)])/0.1)*0.1;
plot([0 maxLim],[0 maxLim],'Color',[0.8 0.8 0.8],'LineWidth',2);
scatter(tempBefFlat(:),tempProbeFlat(:),round(tempSize(~nanFlag)),tempColor(~nanFlag,:),'filled');
xlim([0 maxLim]); ylim([0 maxLim])
xlabel('Bias Reinf BEF'); ylabel('Mean Bias Probe'); title('Bias (Abs) BEF')

subplot(1,3,3); hold on
maxLim = round(max([tempAftFlat(:);tempAftFlat(:)])/0.1)*0.1;
plot([0 maxLim],[0 maxLim],'Color',[0.8 0.8 0.8],'LineWidth',2);
scatter(tempAftFlat(:),tempProbeFlat(:),round(tempSize(~nanFlag)),tempColor(~nanFlag,:),'filled');
xlim([0 maxLim]); ylim([0 maxLim])
xlabel('Bias Reinf AFT'); ylabel('Mean Bias Probe'); title('Bias (Abs) AFT')


%% ALL FUNCTIONS
function [accuracy,accuracyTop,accuracyMean,acc_L,accMean_L,acc_R,accMean_R,...
    bias,bias_absMean,missBias,probeData,reinfDataBef, reinfDataAft,trialNum,probeTrialNum,probeDayNum,dayLen] = getDataByAnimal(mouse, selectProtocol,plotFlag)
global dataPath sep;
load([dataPath sep mouse '.mat']);
    
trialData = fn_readStructByFieldKey(trialData,'trainingType',selectProtocol);

learningCurveBin = 50;

% convert stimulus type from 1 or 2 to 1 or -1 for plotting
stimulus = trialData.stimulus;
stimulus = fn_cell2mat(cellfun(@(x)(smoothdata(-x*2+3,'movmean',...
    learningCurveBin)),stimulus,'UniformOutput',false),1);

responseType = trialData.responseType; 

actionRate = cellfun(@(x)(smoothdata(double(x~=0),'movmean',...
    learningCurveBin)),responseType,'UniformOutput',false);

accuracy = cellfun(@(x)(smoothdata(double(x==1),'movmean',...
    learningCurveBin)),responseType,'UniformOutput',false);

accuracyCellByDay = cellfun(@(x,y)(x./y),accuracy,actionRate,'UniformOutput',false);
accuracyTop = fn_cell2mat(cellfun(@(x)(prctile(x,90)),accuracyCellByDay,'UniformOutput',false),1);
accuracyMean = fn_cell2mat(cellfun(@nanmean,accuracyCellByDay,'UniformOutput',false),1);

[actionRate, ~] = fn_cell2mat(actionRate,1);
[accuracy, cellSize] = fn_cell2mat(accuracy,1);
accuracy = accuracy./actionRate;

[bias, acc_L, acc_R] = cellfun(@getBias,trialData.stimulus,trialData.responseType,'UniformOutput',false);

bias_absMean = fn_removeEmptyEntryCell(bias);
bias_absMean = cellfun(@(x)(nanmean(abs(x))), bias_absMean,'UniformOutput',false);
bias_absMean = fn_cell2mat(bias_absMean,1);

accMean_L = fn_removeEmptyEntryCell(acc_L); accMean_R = fn_removeEmptyEntryCell(acc_R);
accMean_L = fn_cell2mat(cellfun(@nanmean, accMean_L,'UniformOutput',false),1);
accMean_R = fn_cell2mat(cellfun(@nanmean, accMean_R,'UniformOutput',false),1);

bias = fn_cell2mat(bias,1); acc_L = fn_cell2mat(acc_L,1); acc_R = fn_cell2mat(acc_R,1);



missBias = cellfun(@getMissBias,trialData.stimulus,trialData.responseType,'UniformOutput',false);
missBias = fn_cell2mat(missBias,1);
%--------PROBE ORGANIZE BY DAY---------
[probeData, trialNum,probeTrialNum,trialNumNoMiss,probeTrialNumNoMiss,reinfDataBef, reinfDataAft]...
    = cellfun(@getProbe,trialData.stimulus,trialData.responseType,trialData.context,'UniformOutput',false);

probeData = fn_cell2mat(probeData,1);
reinfDataBef = fn_cell2mat(reinfDataBef,1);
reinfDataAft = fn_cell2mat(reinfDataAft,1);

trialNum = fn_cell2mat(trialNum,1);
probeTrialNum = fn_cell2mat(probeTrialNum,1);
probeDayNum = find(~isnan(probeTrialNum));

cumsumTrials = [0;cumsum(trialNum)];
probeTrialNum = probeTrialNum + cumsumTrials(1:end-1);
probeTrialNum(isnan(probeTrialNum)) = [];

dayLen = cellSize(1,:); dayLen(dayLen==0) = [];

if plotFlag
    plotAnimalByDay([accuracy actionRate stimulus bias missBias],probeData,probeTrialNum,dayLen,...
        {'Accuracy','ActionRate','StimulusFreq','ActionBias','MissBias'},mouse)
    
end

end

function plotAnimalByDay(mat,probeData,probeTrialNum,dayLen,ylabels,mouse)

if nargin == 2; ylabels = {[],[]};end
cumsumDayLen = [0 cumsum(dayLen)];

f = figure; set(f,'Units','Normalized','OuterPosition',[0.05,0.58,0.32,0.40])
for i = 1:size(mat,2)
    subplot_tight(size(mat,2),1,i,[0.05 0.05]); hold on;
    if i == 1 || i==2; ylimm = [0 1]; else; ylimm = [-1 1]; end
    if i == 1; plot([1 cumsumDayLen(end)],[0.5 0.5],...
            'Color',[0.8 0.8 0.8],'LineWidth',2);
    elseif i==3 || i==4 || i==5; plot([1 cumsumDayLen(end)],[0 0],...
            'Color',[0.8 0.8 0.8],'LineWidth',2); 
    end
    
    for j = 1:length(dayLen)
        plot([cumsumDayLen(j+1) cumsumDayLen(j+1)], ylimm, 'Color', [0.8 0.8 0.8], 'LineWidth',1.5);
        plot(cumsumDayLen(j)+1:cumsumDayLen(j+1),mat(cumsumDayLen(j)+1:cumsumDayLen(j+1),i),...
            'Color',matlabColors(1,0.8),'LineWidth',2);             
    end   
    
    if i == 1; plot(probeTrialNum,mean(probeData(:,1:2),2),'--o',...
            'Color',[0.2 0.2 0.2],'LineWidth',2);
    elseif i == 2; plot(probeTrialNum,mean(probeData(:,3:4),2),'--o',...
            'Color',[0.2 0.2 0.2],'LineWidth',2);
    elseif i == 4; plot(probeTrialNum,probeData(:,5),'--o',...
            'Color',[0.2 0.2 0.2],'LineWidth',2);
    end      
    
    %xlim([1 2400])
    xlim([1 cumsumDayLen(end)])
    ylim(ylimm);
    ylabel(ylabels{i})
    if i~= size(mat,2); xticks([]); else; xlabel('Trials'); end
    if i==1; title(mouse); end
end


end

function [probeData, trialNum, probeTrialNum, trialNumNoMiss, probeTrialNumNoMiss,...
    reinfDataBef, reinfDataAft] = getProbe(stim, response,ctxt)
    
probeIdx = find(ctxt==3); nProbeRaw = length(probeIdx);

trialNumNoMissIdx = cumsum((response~=0));

probeData = []; reinfDataBef = []; reinfDataAft = [];
goodProbeFlag = true;
% Two probe sessions 
if ~isempty(probeIdx)
    if mod(nProbeRaw,10) == 2 && nProbeRaw==22
        probeIdx([1,12]) = [];
        reinfBefIdx = probeIdx-10; reinfAftIdx = probeIdx+10;
    elseif mod(nProbeRaw,10) == 1 && (nProbeRaw==11 || nProbeRaw==21)
        probeIdx(1) = [];
        if nProbeRaw == 11
            reinfBefIdx = probeIdx-10; reinfAftIdx = probeIdx+10;
        elseif nProbeRaw == 21
            reinfBefIdx = probeIdx-20; reinfAftIdx = probeIdx+20;
        end
        
    elseif nProbeRaw==20
        if (nProbeRaw(end)- nProbeRaw(1))~=19
            reinfBefIdx = probeIdx-10; reinfAftIdx = probeIdx+10;
        else 
            reinfBefIdx = probeIdx-20; reinfAftIdx = probeIdx+20;
        end   
    elseif nProbeRaw==10
        reinfBefIdx = probeIdx-10; reinfAftIdx = probeIdx+10;
    else
        disp(['WARNING -- Exception in Probe! nProbe = ' int2str(nProbeRaw)]);
        probeTrialNum = nan; probeTrialNumNoMiss = nan;
        goodProbeFlag = false;
    end
    if goodProbeFlag
        probeData = getIdxAccBias(stim, response, probeIdx);
        reinfDataBef = getIdxAccBias(stim, response, reinfBefIdx);
        reinfDataAft = getIdxAccBias(stim, response, reinfAftIdx);

        probeTrialNum = round(mean(probeIdx));

        probeTrialNumNoMiss = round(mean(trialNumNoMissIdx(probeIdx)));
    end
else
    probeTrialNum = nan; probeTrialNumNoMiss = nan;
end
if isempty(ctxt); trialNum = []; trialNumNoMiss = []; %probeTrialNum = []; probeTrialNumNoMiss = []; 
else; trialNum = length(ctxt); trialNumNoMiss = max(trialNumNoMissIdx);
end 

end

function idxData = getIdxAccBias(stim, resp, idx)

if max(idx) <= length(stim)
    idxStim = stim(idx); idxResponse = resp(idx);
    tempResp = idxResponse(idxStim==1);

    idxData(1) = sum(tempResp==1) / sum(tempResp~=0); % Accuracy L
    idxData(3) = sum(tempResp~=0) / length(tempResp); % Action rate L

    tempResp = idxResponse(idxStim==2);
    idxData(2) = sum(tempResp==1) / sum(tempResp~=0); % Accuracy R
    idxData(4) = sum(tempResp~=0) / length(tempResp); % Action rate R

    idxData(5) = (idxData(1) - idxData(2)); %/ (probeData(1) + probeData(2));
else
    idxData = nan(1,5);
end

end

function probeData = getProbePerformance(stim, response,probeStart,probeEnd)
probeData = zeros(length(probeStart),5);
for i = 1:length(probeStart)
    tempStim = stim(probeStart(i):probeEnd(i));
    tempResp = response(probeStart(i):probeEnd(i));

    tempRespS1 = tempResp(tempStim==1);
    probeData(i,1) = sum(tempRespS1==1) / sum(tempRespS1~=0); % Accuracy L, miss corrected
    probeData(i,3) = sum(tempRespS1~=0) / length(tempRespS1); % Action rate L

    tempRespS2 = tempResp(tempStim==2);
    probeData(i,2) = sum(tempRespS2==1) / sum(tempRespS2~=0); % Accuracy R, miss corrected
    probeData(i,4) = sum(tempRespS2~=0) / length(tempRespS2); % Action rate R   
end
probeData(i,5) = (probeData(i,1) - probeData(i,2)) ./ (probeData(i,1) + probeData(i,2));
end

function plotAllAnimalsReinfProbe(accuracyCell,probeCell,probeTrialNumCell)
trialLim = 2400;
trialBin = 0:300:trialLim;
probeAllAni = nan(length(probeCell),length(trialBin)-1);

for i= 1:length(probeCell)
    trialNum = probeTrialNumCell{i}; 
    trialLimFlag = trialNum>trialLim; trialNum(trialLimFlag) = [];
    tempPerf = probeCell{i}(~trialLimFlag,:);
    tempPerf = mean(tempPerf(:,1:2),2);
    
    binIdx = sum(trialNum>trialBin,2);
    probeAllAni(i,binIdx) = tempPerf;
end

meanProbePerf = nanmean(probeAllAni,1); 
nData = sum(~isnan(probeAllAni),1);
semProbePerf = nanstd(probeAllAni,0,1)./sqrt(nData);

nanflag = isnan(meanProbePerf);
trialBinPlot = (trialBin(1:end-1) + trialBin(2:end))/2;

% REINF DATA
accuracyCell = cellfun(@(x)(x(1:trialLim)),accuracyCell,'UniformOutput',false);
accuracyCell = fn_cell2mat(accuracyCell,2);

% PLOT 1
figure;
hold on;
f_errorbar = fn_plotFillErrorbar(trialBinPlot(~nanflag),meanProbePerf(~nanflag),semProbePerf(~nanflag),...
    matlabColors(2),'faceAlpha',0.2,'LineStyle','none');
f_errorbar = fn_plotFillErrorbar(1:size(accuracyCell,1),(mean(accuracyCell,2))',...
    (std(accuracyCell,0,2)./sqrt(size(accuracyCell,2)))',...
    matlabColors(1),'faceAlpha',0.2,'LineStyle','none');
% PLOT REINF
plot(mean(accuracyCell,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
% PLOT PROBE AGAIN
plot(trialBinPlot(~nanflag),meanProbePerf(~nanflag),'--o','Color',matlabColors(2,0.9),'LineWidth',2)

% PLOT 2
figure;
hold on;
% PLOT REINF
for j = 1:size(accuracyCell,2)
    plot(accuracyCell(:,j),'Color', matlabColors(1,0.2), 'LineWidth',1);
end
for j = 1:size(probeAllAni,1)
    plot(trialBinPlot,probeAllAni(j,:),'--o','Color', matlabColors(2,0.4), 'LineWidth',1)
end
plot(mean(accuracyCell,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
% PLOT PROBE AGAIN
plot(trialBinPlot(~nanflag),meanProbePerf(~nanflag),'--o','Color',matlabColors(2,0.9),'LineWidth',2)

end

function plotAllAnimalsReinfProbeBias(biasCell,probeCell,probeTrialNumCell)
trialLim = 2400;
trialBin = 0:300:trialLim;
probeAllAni = nan(length(probeCell),length(trialBin)-1);

for i= 1:length(probeCell)
    trialNum = probeTrialNumCell{i}; 
    trialLimFlag = trialNum>trialLim; trialNum(trialLimFlag) = [];
    tempPerf = probeCell{i}(~trialLimFlag,:);
    tempPerf = abs(tempPerf(:,5));
    
    binIdx = sum(trialNum>trialBin,2);
    probeAllAni(i,binIdx) = tempPerf;
end

meanProbePerf = nanmean(probeAllAni,1); 
nData = sum(~isnan(probeAllAni),1);
semProbePerf = nanstd(probeAllAni,0,1)./sqrt(nData);

nanflag = isnan(meanProbePerf);
trialBinPlot = (trialBin(1:end-1) + trialBin(2:end))/2;

% REINF DATA
biasCell = cellfun(@(x)(x(1:trialLim)),biasCell,'UniformOutput',false);
biasCell = fn_cell2mat(biasCell,2);

% PLOT 1
figure;
hold on;
f_errorbar = fn_plotFillErrorbar(trialBinPlot(~nanflag),meanProbePerf(~nanflag),semProbePerf(~nanflag),...
    matlabColors(2),'faceAlpha',0.2,'LineStyle','none');
tempMean =(nanmean(biasCell,2))'; tempSEM = (nanstd(biasCell,0,2)./sqrt(size(biasCell,2)))';
tempX = 1:size(biasCell,1); tempNanFlag = ~isnan(tempMean);
f_errorbar = fn_plotFillErrorbar(tempX(tempNanFlag),tempMean(tempNanFlag),...
    tempSEM(tempNanFlag),matlabColors(1),'faceAlpha',0.2,'LineStyle','none');
% PLOT REINF
plot(smoothdata(nanmean(biasCell,2),'movmean',50),'Color', matlabColors(1,0.9), 'LineWidth',2);
% PLOT PROBE AGAIN
plot(trialBinPlot(~nanflag),meanProbePerf(~nanflag),'--o','Color',matlabColors(2,0.9),'LineWidth',2)

% PLOT 2
figure;
hold on;
% PLOT REINF
for j = 1:size(biasCell,2)
    plot(biasCell(:,j),'Color', matlabColors(1,0.2), 'LineWidth',1);
end
for j = 1:size(probeAllAni,1)
    plot(trialBinPlot,probeAllAni(j,:),'--o','Color', matlabColors(2,0.4), 'LineWidth',1)
end
plot(smoothdata(nanmean(biasCell,2),'movmean',50),'Color', matlabColors(1,0.9), 'LineWidth',2);
% PLOT PROBE AGAIN
plot(trialBinPlot(~nanflag),meanProbePerf(~nanflag),'--o','Color',matlabColors(2,0.9),'LineWidth',2)

end

function [bias,acc_L,acc_R] = getBias(stimulus,responeType)
    biasBin = 100;
    correct = (responeType == 1);
    miss = (responeType == 0);
    bias = nan(size(correct)); acc_L = nan(size(correct)); acc_R = nan(size(correct));
    for j = biasBin:length(correct)
        tempCorrect = correct((j-biasBin+1):j);
        tempMiss = miss((j-biasBin+1):j);
        tempStim = stimulus((j-biasBin+1):j);
        acc_s1 = sum(tempCorrect & tempStim==1) / sum(tempStim==1 & tempMiss~=1);
        acc_s2 = sum(tempCorrect & tempStim==2) / sum(tempStim==2 & tempMiss~=1);
        bias(j - biasBin/2) = (acc_s1 - acc_s2); % / (acc_s1 + acc_s2);
        acc_L(j - biasBin/2) = acc_s1; acc_R(j - biasBin/2) = acc_s2;
    end
end

function bias = getMissBias(stimulus,responeType)
    biasBin = 50;
    miss = (responeType == 0);
    bias = nan(size(miss));
    for j = biasBin:length(miss)
        tempMiss = miss((j-biasBin+1):j);
        tempStim = stimulus((j-biasBin+1):j);
        acc_s1 = sum(tempMiss & tempStim==1) / sum(tempStim==1);
        acc_s2 = sum(tempMiss & tempStim==2) / sum(tempStim==2);
        if sum(tempMiss) < biasBin*0.2
            bias(j - biasBin/2) = nan;
        else
            bias(j - biasBin/2) = (acc_s1 - acc_s2) / (acc_s1 + acc_s2);
        end
    end
end
