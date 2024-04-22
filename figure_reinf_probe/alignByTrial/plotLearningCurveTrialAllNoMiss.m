function plotLearningCurveTrialAllNoMiss(mice,selectProtocol)

if nargin == 1; selectProtocol = {'puretone'}; end
global sep;

rootPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\';
loadPath = [rootPath sep 'trialData\'];
figPath = [rootPath sep 'figure\' sep 'learningCurveTrialAll' sep]; mkdir(figPath);
trialData = cellfun(@(x)fn_selectProtocol(loadPath, x, selectProtocol),mice,'UniformOutput',false);

%[~, accuracy,~,~,~,~,~,~,~,~,~,~,probeData,~, ~,~,probeTrialNum,~,~]
outmat = cellfun(@(x)alignDays(x),trialData,'UniformOutput',false);

accuracy = cellfun(@(x)(x.accuracy),outmat,'UniformOutput',false);
bias = cellfun(@(x)(x.bias),outmat,'UniformOutput',false);
probeData = cellfun(@(x)(x.probeData),outmat,'UniformOutput',false);
probeTrialNum = cellfun(@(x)(x.probeTrialNum),outmat,'UniformOutput',false);
reactionTime = cellfun(@(x)(x.reactionTime),outmat,'UniformOutput',false);
lowActionDayFlag = cellfun(@(x)(x.lowActionDayFlag),outmat,'UniformOutput',false);
wheelPreSound = cellfun(@(x)(x.wheelPreSound),outmat,'UniformOutput',false);
wheelSoundOn = cellfun(@(x)(x.wheelSoundOn),outmat,'UniformOutput',false);
reinfDataBef = cellfun(@(x)(x.reinfDataBef),outmat,'UniformOutput',false);
reinfDataAft = cellfun(@(x)(x.reinfDataAft),outmat,'UniformOutput',false);
probeDayNum = cellfun(@(x)(x.probeDayNum),outmat,'UniformOutput',false);

trialLim = 3000;
trialBin = 0:250:trialLim;
probeAllAni = nan(length(probeData),length(trialBin)-1);
probeAllAniL = nan(length(probeData),length(trialBin)-1);
probeAllAniR = nan(length(probeData),length(trialBin)-1);
probeAllAniBias = nan(length(probeData),length(trialBin)-1);
for i= 1:length(probeData)
    trialNum = probeTrialNum{i}; 
    trialLimFlag = trialNum>trialLim; trialNum(trialLimFlag) = [];
    tempData = probeData{i}(~trialLimFlag,:);
    tempPerf = mean(tempData(:,1:2),2);
    tempBias = abs(tempData(:,5));
    
    binIdx = sum(trialNum>trialBin,2);
    probeAllAni(i,binIdx) = tempPerf;
    probeAllAniL(i,binIdx) = tempData(:,1);
    probeAllAniR(i,binIdx) = tempData(:,2);
    probeAllAniBias(i,binIdx) = tempBias;
end


nData = sum(~isnan(probeAllAni),1);
meanProbePerf = nanmean(probeAllAni,1); 
semProbePerf = nanstd(probeAllAni,0,1)./sqrt(nData);
meanProbeBias = nanmean(probeAllAniBias,1); 
semProbeBias = nanstd(probeAllAniBias,0,1)./sqrt(nData);

nanflag = isnan(meanProbePerf);
trialBinPlot = (trialBin(1:end-1) + trialBin(2:end))/2;

% REINF DATA
accuracy = cellfun(@(x)fn_attachNan(x,trialLim,1), accuracy, 'UniformOutput',false);
accuracy = cellfun(@(x)(x(1:trialLim)),accuracy,'UniformOutput',false);
accuracy = fn_cell2mat(accuracy,2);

%--------------------------------PLOT 1-----------------------------------
fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
f_errorbar = fn_plotFillErrorbar(trialBinPlot(~nanflag),meanProbePerf(~nanflag),semProbePerf(~nanflag),...
    matlabColors(2),'faceAlpha',0.2,'LineStyle','none');
f_errorbar = fn_plotFillErrorbar(1:size(accuracy,1),(nanmean(accuracy,2))',...
    (nanstd(accuracy,0,2)./sqrt(size(accuracy,2)))',...
    matlabColors(1),'faceAlpha',0.2,'LineStyle','none');
% PLOT REINF
plot(nanmean(accuracy,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
% PLOT PROBE AGAIN
plot(trialBinPlot(~nanflag),meanProbePerf(~nanflag),'--o','Color',matlabColors(2,0.9),'LineWidth',2)
plot([1 trialLim],[0.5 0.5],'Color',[0.6 0.6 0.6],'LineWidth',2);
xlim([1 trialLim]); ylabel('Accuracy'); xlabel('Trials');

%--------------------------------PLOT 2-----------------------------------
fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
% PLOT REINF
for j = 1:size(accuracy,2)
    plot(accuracy(:,j),'Color', matlabColors(1,0.2), 'LineWidth',1);
end
for j = 1:size(probeAllAni,1)
    plot(trialBinPlot,probeAllAni(j,:),'--o','Color', matlabColors(2,0.4), 'LineWidth',1)
end
plot(nanmean(accuracy,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
% PLOT PROBE AGAIN
plot(trialBinPlot(~nanflag),meanProbePerf(~nanflag),'--o','Color',matlabColors(2,0.9),'LineWidth',2)
plot([1 trialLim],[0.5 0.5],'Color',[0.6 0.6 0.6],'LineWidth',2);
xlim([1 trialLim]); ylabel('Accuracy'); xlabel('Trials');

%--------------------------------PLOT 3-----------------------------------
fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
% PLOT REINF
for j = 1:size(accuracy,2)
    plot(accuracy(:,j),'Color', matlabColors(1,0.2), 'LineWidth',1);
end
plot(nanmean(accuracy,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
plot([1 trialLim],[0.5 0.5],'Color',[0.6 0.6 0.6],'LineWidth',2);
xlim([1 trialLim]); ylabel('Accuracy'); xlabel('Trials');

%-----------------------REINF BEF AFT WITH PROBE--------------------------
nDay = 8;
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

nDay = 7;
probeBiasPlot = abs(probeDayAlign(1:nDay,:,5));
reinfBefBiasPlot = abs(reinfBefDayAlign(1:nDay,:,5));
reinfAftBiasPlot = abs(reinfAftDayAlign(1:nDay,:,5));

figure; subplot(1,2,1); hold on;
tempBefFlat = reshape(reinfBefDayPlot(1:nDay,:),1,[]);
tempAftFlat = reshape(reinfAftDayPlot(1:nDay,:),1,[]);
tempProbeFlat = reshape(probeDayPlot(1:nDay,:),1,[]);
nanFlag = isnan(tempProbeFlat); 
tempBefFlat(nanFlag) = []; tempAftFlat(nanFlag) = []; tempProbeFlat(nanFlag) = [];
[pBef,hBef] = signrank(tempBefFlat,tempProbeFlat,'tail','left');
[pAft,hAft] = signrank(tempAftFlat,tempProbeFlat,'tail','left');
bar([nanmean(tempBefFlat) nanmean(tempProbeFlat) nanmean(tempAftFlat)] ,'EdgeColor','None','FaceColor',matlabColors(2,0.9)); hold on;
for i = 1:length(tempBefFlat)
    f = plot([1 2 3],[tempBefFlat(i) tempProbeFlat(i) tempAftFlat(i)],'Color',[0.6 0.6 0.6],'Marker','.','MarkerSize',15,...
        'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',0.5);
end
legend(f,['pBef = ' num2str(pBef,'%.2e') newline 'pAft = ' num2str(pAft,'%.2e')],'Location','Best')
ylabel('Accuracy'); xticks([1 2 3]); xticklabels({'Bef','Probe','Aft'}); xlim([0 4]);

subplot(1,2,2); hold on;
tempBefFlat = reshape(reinfBefBiasPlot(1:nDay,:),1,[]);
tempAftFlat = reshape(reinfAftBiasPlot(1:nDay,:),1,[]);
tempProbeFlat = reshape(probeBiasPlot(1:nDay,:),1,[]);
nanFlag = isnan(tempProbeFlat); 
tempBefFlat(nanFlag) = []; tempAftFlat(nanFlag) = []; tempProbeFlat(nanFlag) = [];
[pBef,hBef] = signrank(tempBefFlat,tempProbeFlat,'tail','right');
[pAft,hAft] = signrank(tempAftFlat,tempProbeFlat,'tail','right');
bar([nanmean(tempBefFlat) nanmean(tempProbeFlat) nanmean(tempAftFlat)] ,'EdgeColor','None','FaceColor',matlabColors(2,0.9)); hold on;
for i = 1:length(tempBefFlat)
    f = plot([1 2 3],[tempBefFlat(i) tempProbeFlat(i) tempAftFlat(i)],'Color',[0.6 0.6 0.6],'Marker','.','MarkerSize',15,...
        'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',0.5);
end
legend(f,['pBef = ' num2str(pBef,'%.2e') newline 'pAft = ' num2str(pAft,'%.2e')],'Location','Best')
ylabel('Bias(Abs)'); xticks([1 2 3]); xticklabels({'Bef','Probe','Aft'}); xlim([0 4]);

%------------------------------REACTION TIME-------------------------------
reactionTime = cellfun(@(x,y)(fn_fillMat(x,logical(y))),reactionTime,lowActionDayFlag,'UniformOutput',false);
reactionTime = cellfun(@(x)fn_attachNan(x,trialLim,1), reactionTime, 'UniformOutput',false);
reactionTime = cellfun(@(x)(x(1:trialLim)),reactionTime,'UniformOutput',false);
reactionTime = fn_cell2mat(reactionTime,2);
reactionTime(reactionTime > 2.6) = nan; reactionTime(reactionTime < 0) = nan;
reactionTime = smoothdata(reactionTime,1,'movmean',50,'includenan');

fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
for j = 1:size(accuracy,2)
    plot(reactionTime(:,j),'Color', matlabColors(1,0.2), 'LineWidth',1);
end
plot(nanmean(reactionTime,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
ylim([0.3 1.5]); xlim([1 trialLim]);
ylabel('Reaction Time (s)'); xlabel('Trials');

%------------------------------Bias-------------------------------
bias = cellfun(@(x)fn_attachNan(x,trialLim,1), bias, 'UniformOutput',false);
bias = cellfun(@(x)(x(1:trialLim)),bias,'UniformOutput',false);
bias = fn_cell2mat(bias,2);
bias = abs(bias);
% REINF ONLY
fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
for j = 1:size(bias,2)
    plot(bias(:,j),'Color', matlabColors(1,0.2), 'LineWidth',1);
end
plot(nanmean(bias,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
xlim([1 trialLim]); ylim([0 1]); ylabel('Bias'); xlabel('Trials');
% REINF AND PROBE
fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
plot(nanmean(bias,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
f_errorbar = fn_plotFillErrorbar(1:size(bias,1),(nanmean(bias,2))',...
    (nanstd(bias,0,2)./sqrt(size(bias,2)))',...
    matlabColors(1),'faceAlpha',0.2,'LineStyle','none');
f_errorbar = fn_plotFillErrorbar(trialBinPlot(~nanflag),meanProbeBias(~nanflag),semProbeBias(~nanflag),...
    matlabColors(2),'faceAlpha',0.2,'LineStyle','none');
plot(trialBinPlot(~nanflag),meanProbeBias(~nanflag),'--o','Color',matlabColors(2,0.9),'LineWidth',2)
xlim([1 trialLim]); ylim([0 1]); ylabel('Bias'); xlabel('Trials');

%---------------------PLOT ACCURACY AND BIAS TOGETHER----------------------
fn_figureSmartDim('hSize',0.3,'widthHeightRatio',1) ; hold on; 
avgAccSmooth = smoothdata(nanmean(accuracy,2),'movmean',200);
avgBiasSmooth = smoothdata(nanmean(bias,2),'movmean',200);
for i = 1:size(accuracy,2)
    accSmooth = smoothdata(accuracy(:,i),'movmean',400); biasSmooth = smoothdata(bias(:,i),'movmean',400);
    plot(accSmooth,biasSmooth,'Color',matlabColors(i,0.3));
    plot(accSmooth(1),biasSmooth(1),'.','MarkerSize',15,'Color',matlabColors(i))
    plot(accSmooth(end),biasSmooth(end),'*','MarkerSize',8,'Color',matlabColors(i))
end
plot(avgAccSmooth,avgBiasSmooth,'Color',[0 0 0],'LineWidth',2);
xlim([0.3 0.95]); xticks([0.3 0.5 0.7 0.9])
xlabel('Accuracy'); ylabel('Action Bias')

%------------------------------ITI MOVEMENTS-------------------------------
wheelPreSound = cellfun(@(x)fn_attachNan(x,trialLim,1), wheelPreSound, 'UniformOutput',false);
wheelPreSound = cellfun(@(x)(x(1:trialLim)),wheelPreSound,'UniformOutput',false);
wheelPreSound = fn_cell2mat(wheelPreSound,2);
wheelPreSound = smoothdata(wheelPreSound,1,'movmean',50,'includenan');
%wheelPreSound = wheelPreSound./repmat(nanmean(wheelPreSound(1:500,:),1),trialLim,1); %normalize to the first day of each animal
fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
for j = 1:size(accuracy,2)
    plot(wheelPreSound(:,j),'Color', matlabColors(1,0.2), 'LineWidth',1);
end
plot(nanmean(wheelPreSound,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
xlim([1 trialLim]);
ylabel('Reaction Time (s)'); xlabel('Trials');

%-------------------ITI MOVEMENTS vs. Sound-on Movement--------------------
wheelSoundOn = cellfun(@(x)fn_attachNan(x,trialLim,1), wheelSoundOn, 'UniformOutput',false);
wheelSoundOn = cellfun(@(x)(x(1:trialLim)),wheelSoundOn,'UniformOutput',false);
wheelSoundOn = fn_cell2mat(wheelSoundOn,2);

fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
temp = [nanmean(wheelPreSound,1);  nanmean(wheelSoundOn,1)] * 0.35;
bar(nanmean(temp,2),'EdgeColor',[1 1 1],'FaceColor',matlabColors(2,0.9)); hold on;
for i = 1:size(temp,1)
    scatter(i*ones(size(temp,2),1) ,temp(i,:),10,[0.6 0.6 0.6],'filled');
end
plot(temp,'Color',[0.6 0.6 0.6],'LineWidth',0.8)
ylabel('Degree of Movement'); xticks([1 2]);xticklabels({'ITI',...
    sprintf('%s\\newline%s\n', 'Response','Window')});
xtickangle(0);

%---------------------------PLOT MODEL TOGETHER----------------------------
weightCell = readModelPerformance(mice);
if ~isempty(weightCell)
    modelPerf = cellfun(@(x)fn_attachNan(x,trialLim,2), weightCell, 'UniformOutput',false);
    modelPerf = cellfun(@(x)(x(1:trialLim)),modelPerf,'UniformOutput',false);
    modelPerf = fn_cell2mat(modelPerf,1)';

    fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
    % PLOT REINF
    plot(nanmean(accuracy,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
    f_errorbar = fn_plotFillErrorbar(1:size(accuracy,1),(nanmean(accuracy,2))',...
        (nanstd(accuracy,0,2)./sqrt(size(accuracy,2)))',...
        matlabColors(1),'faceAlpha',0.2,'LineStyle','none');
    % PLOT PROBE 
    plot(trialBinPlot(~nanflag),meanProbePerf(~nanflag),'--o','Color',matlabColors(2,0.9),'LineWidth',2)
    f_errorbar = fn_plotFillErrorbar(trialBinPlot(~nanflag),meanProbePerf(~nanflag),semProbePerf(~nanflag),...
        matlabColors(2),'faceAlpha',0.2,'LineStyle','none');

    % PLOT MODEL
    plot(nanmean(modelPerf,2),'Color', matlabColors(7,0.9), 'LineWidth',2);
    f_errorbar = fn_plotFillErrorbar(1:size(modelPerf,1),(nanmean(modelPerf,2))',...
        (nanstd(modelPerf,0,2)./sqrt(size(modelPerf,2)))',...
        matlabColors(7),'faceAlpha',0.2,'LineStyle','none');
    plot([1 trialLim],[0.5 0.5],'Color',[0.6 0.6 0.6],'LineWidth',2);
    xlim([1 trialLim]); ylabel('Accuracy'); xlabel('Trials');

    %a = cellfun(@(x)(find(mean(x(:,1:2),2)>0.75)),probeData,'UniformOutput',false);
    %b = cellfun(@(x)(x(2)),a,'UniformOutput',false);
    %probeThreTrialNum = cellfun(@(x,y)(x(y)),probeTrialNum,b);
    %for i = 1:size(bias,2)
    %   bias_befProbe(i) = nanmean(bias(1:probeThreTrialNum,i));  
    %
    %end
end

%----------------------------SOME TESTING CODE-----------------------------
% meanAcc = smoothdata(accuracy,'movmean',500);
% meanBias = smoothdata(bias,'movmean',500);
% for i = 1:7; firstHit70Temp = find(meanAcc(:,i)>0.7);
%     if ~isempty(firstHit70Temp); firstHit70(i) = firstHit70Temp(1);
%     else firstHit70(i) = inf;
%     end
% end
% for i = 1:7; firstBias03Temp = find(meanBias(400:end,i)<0.3);
%     if ~isempty(firstBias03Temp); firstBias03(i) = firstBias03Temp(1)+400;
%     else firstBias03(i) = inf;
%     end
% end

%figure;
%scatter(max(firstHit70,firstBias03),nanmean(bias(1:700,:),1))

end