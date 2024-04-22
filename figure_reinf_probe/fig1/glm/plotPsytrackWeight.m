clear;
global loadPath;
loadPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\psyTrackFit\';
%modelComparison = {'S','SA','SRp','SRn','SRpRn','SB','SBA','SBRp','SBRn','SBRpRn'};

%% -------------------LOAD WEIGHTS------------------------
nPrev = 1;
mouse ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
modelComparison = 'SBARpRn';
learningLim = [];

sW = []; bW = []; ahW = []; rnW = []; rpW = [];
for i = 1:length(mouse)
    load([loadPath filesep mouse{i} 'psytrack_' modelComparison '_nPrev' int2str(nPrev) '.mat']);
    sW{i} = abs(wMode(5,:)); bW{i} = abs(wMode(4,:)); ahW{i} = abs(wMode(1,:)); rnW{i} = abs(wMode(2,:)); rpW{i} = abs(wMode(3,:));
end
sW = fn_cell2matFillNan(sW); bW = fn_cell2matFillNan(bW); ahW = fn_cell2matFillNan(ahW); rnW = fn_cell2matFillNan(rnW);rpW = fn_cell2matFillNan(rpW);

if ~isnan(learningLim)
    sW = sW(:,1:learningLim); bW = bW(:,1:learningLim); ahW = ahW(:,1:learningLim); rnW = rnW(:,1:learningLim); rpW = rpW(:,1:learningLim);
end
try  sW(12,1738:end) = nan; bW(12,1738:end) = nan; ahW(12,1738:end) = nan; rnW(12,1738:end) = nan; rpW(12,1738:end) = nan;
catch; end
%% -------------------VISUALIZE WEIGHT GAUSSIAN GRAPH------------------------

earlyLearning = 1:100; lateLearning = 1801:2000;
figure; 
subplot(1,2,1); 
fn_plotBarPaired({nanmean(sW(:,earlyLearning),2), nanmean(sW(:,lateLearning),2)});
subplot(1,2,2);
earlyLearning = 601:800; lateLearning = 1801:2000;
fn_plotBarPaired({nanmean(bW(:,earlyLearning),2), nanmean(bW(:,lateLearning),2)});

tempPlot = {nanmean(sW,2),nanmean(bW,2),nanmean(ahW,2),nanmean(rpW,2),nanmean(rnW,2)};
figure; hold on;
plot([0 6],[0 0],'LineWidth',2,'Color',[0.8 0.8 0.8])
violin(tempPlot,'facecolor',matlabColors([1:4 7]),'edgecolor',[0.4 0.4 0.4],'mc',[0.4 0.4 0.4],'medc',[],'facealpha',0.5);

for i = 1:length(tempPlot)
    scatter(ones(1,length(tempPlot{i}))*i,tempPlot{i},8,[0.4 0.4 0.4],'filled');
end
xlim([0.5 5.5]); xticks(1:5); xticklabels({'s';'b';'ah';'rph';'rnh'})

% do repeated anova
%wMat = fn_cell2mat(tempPlot,2); wMat = array2table(wMat,'VariableNames',{'V1','V2','V3','V4','V5'});
%within = table({'V1';'V2';'V3';'V4';'V5'}, 'VariableNames', {'weight'});
% fit the repeated measures model
%rm = fitrm(wMat,'V1-V5~1','WithinDesign',within);
%[ranovatblb] = ranova(rm);
%Mrm2 = multcompare(rm,'HRs','By','ObstaclePos','ComparisonType','bonferroni');

%% -------------------VISUALIZE WEIGHT TOGETHER COSYNE------------------------
sWplot = sW; bWplot = bW; ahWplot = ahW;  rnWplot = rnW;  rpWplot = rpW; 
plotLim = 3000;
for i = 1:length(mouse)
    firstNan = find(isnan(sW(i,:)),1);
    if i ~= 4 && i ~= 6 && i ~= 12 && ~isempty(firstNan)
        sWplot(i,firstNan:end) = sWplot(i,firstNan-1); bWplot(i,firstNan:end) = bWplot(i,firstNan-1);
        ahWplot(i,firstNan:end) = ahWplot(i,firstNan-1); rnWplot(i,firstNan:end) = rnWplot(i,firstNan-1);
        rpWplot(i,firstNan:end) = rpWplot(i,firstNan-1);
    elseif isempty(firstNan); firstNan = size(sW,2);
    end
    %disp([mouse{i} ' trial = ' int2str(firstNan) ' stim = ' num2str(sW(i,firstNan),'%.3f') ', bias = ' num2str(bW(i,firstNan),'%.3f')])
end
sWplot = sWplot(:,1:plotLim);bWplot = bWplot(:,1:plotLim);ahWplot = ahWplot(:,1:plotLim);
rnWplot = rnWplot(:,1:plotLim);rpWplot = rpWplot(:,1:plotLim);

figure; hold on;  plotLim = 3000; 
f_errorbar = fn_plotFillErrorbar(1:size(sWplot,2),nanmean(sWplot,1),nanstd(sWplot,0,1)./sqrt(sum(~isnan(sWplot),1)),...
    matlabColors(1),'faceAlpha',0.2,'LineStyle','none');
plot(nanmean(sWplot,1),'Color', matlabColors(1,0.9), 'LineWidth',2);

f_errorbar = fn_plotFillErrorbar(1:size(bWplot,2),nanmean(bWplot,1),nanstd(bWplot,0,1)./sqrt(sum(~isnan(bWplot),1)),...
    matlabColors(2),'faceAlpha',0.2,'LineStyle','none');
plot(nanmean(bWplot,1),'Color', matlabColors(2,0.9), 'LineWidth',2);

f_errorbar = fn_plotFillErrorbar(1:size(ahWplot,2),nanmean(ahWplot,1),nanstd(ahWplot,0,1)./sqrt(sum(~isnan(ahWplot),1)),...
    matlabColors(3),'faceAlpha',0.2,'LineStyle','none');
plot(nanmean(ahWplot,1),'Color', matlabColors(3,0.9), 'LineWidth',2);

f_errorbar = fn_plotFillErrorbar(1:size(rnWplot,2),nanmean(rnWplot,1),nanstd(rnWplot,0,1)./sqrt(sum(~isnan(rnWplot),1)),...
    matlabColors(4),'faceAlpha',0.2,'LineStyle','none');
plot(nanmean(rnWplot,1),'Color', matlabColors(4,0.9), 'LineWidth',2);

f_errorbar = fn_plotFillErrorbar(1:size(rpWplot,2),nanmean(rpWplot,1),nanstd(rpWplot,0,1)./sqrt(sum(~isnan(rpWplot),1)),...
    matlabColors(7),'faceAlpha',0.2,'LineStyle','none');
plot(nanmean(rpWplot,1),'Color', matlabColors(7,0.9), 'LineWidth',2);

ylim([-0.2 2.5]); yticks(0:0.5:2.5);
xlim([0 plotLim]);xticks(0:1000:plotLim)
disp('none')


%% CORRELATION BETWEEN weights
figure; subplot(1,2,1);hold on; scatter(tempPlot{2},(tempPlot{5}+tempPlot{4}+tempPlot{3})/3,20,'filled')
lm = fitlm(tempPlot{2},(tempPlot{5}+tempPlot{4}+tempPlot{3})/3);
xAxis = 0:0.01:1.5; yAxis = lm.Coefficients{1,1} + lm.Coefficients{2,1} * xAxis;
plot(xAxis,yAxis,'Color',[0.8 0.8 0.8]);
title(['corr = ' num2str(corr(tempPlot{2},(tempPlot{5}+tempPlot{4}+tempPlot{3})/3)) ',p = ' num2str(lm.Coefficients.pValue(2))])


subplot(1,2,2); hold on; scatter(tempPlot{5},tempPlot{4},20,'filled')
lm = fitlm(tempPlot{5},tempPlot{4});
xlim([0 1.2]); ylim([0 0.6]); yticks(0:0.2:0.6); xticks(0:0.4:1.2);
xAxis = 0:0.01:1.2; yAxis = lm.Coefficients{1,1} + lm.Coefficients{2,1} * xAxis;
plot(xAxis,yAxis,'Color',[0.8 0.8 0.8]);
title(['corr = ' num2str(corr(tempPlot{5},tempPlot{4})) ',p = ' num2str(lm.Coefficients.pValue(2))])

%subjects = repmat((1:13)',[5 1]); weights = [ones(13,1); ones(13,1)*2; ones(13,1)*3; ones(13,1)*4; ones(13,1)*5];
%tempPlot = fn_cell2mat(tempPlot,1);
%save('C:\Users\zzhu34\Documents\gitRep\analyze2AFCObj\pythonData\psytrackW.mat','tempPlot','subjects','weights');
%% ANOVA REPEATED EXAMPLE

Biases = rand(18,3*5*2); % subjects, HRs*obstacle pos.*VFs
varNames = cell(3*5*2,1);
for i = 1 : 3*5*2
v = strcat('V',num2str(i));
varNames{i,1} = v;
end
% Create a table storing the respones
tbiases = array2table(Biases, 'VariableNames',varNames);
% Create a table reflecting the within subject factors
HRs = cell(3*5*2,1); % head roll conditions
VFs = cell(3*5*2,1); % Visual feedback conditions
OPs = cell(3*5*2,1); % Obstacle Positions
% Assiging the values to the parameters based on the data sorting
c1 = cell(1,1); c1{1} = 'Y'; c1 = repmat(c1,15,1); VFs(1: 15,1) = c1;
c1 = cell(1,1); c1{1} = 'N'; c1 = repmat(c1,15,1); VFs(16: end,1) = c1;
c1 = cell(1,1); c1{1} = 'HR0'; c1 = repmat(c1,10,1); HRs(1:3:end,1) = c1;
c1 = cell(1,1); c1{1} = 'HRL'; c1 = repmat(c1,10,1); HRs(2:3:end,1) = c1;
c1 = cell(1,1); c1{1} = 'HRR'; c1 = repmat(c1,10,1); HRs(3:3:end,1) = c1;
for i = 1 : 5
o = strcat('O',num2str(i));
c1 = cell(1,1); c1{1} = o; c1 = repmat(c1,3,1); OPs((i-1)*3+1:i*3,1) =c1;
end
OPs(16:end,1) = OPs(1:15,1);
% Create the within table
factorNames = {'HRs','VisualFeedback', 'ObstaclePos'};
within = table(HRs, VFs, OPs, 'VariableNames', factorNames);


%% ----------------MODEL FIT COSYNE-------------------
clear;
global loadPath;
nPrev = 1;
loadPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\psyTrackFit\';

mouse ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
rsq = [];
for i = 1:length(mouse)
    %mouse = 'zz063';
    modelComparison = 'SBARpRn';
    load(['C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\trialData' filesep mouse{i} '_nPrev5.mat']);
    
    load([loadPath filesep mouse{i} 'psytrack_' modelComparison '_nPrev' int2str(nPrev) '.mat']);
    regressor = [actionH(:,1),actionXnegRewardH(:,1),actionXposRewardH(:,1),ones(size(stimulus)),stimulus];
    predChoice = fn_logistic (sum(regressor .* wMode',2));
    predAcc = (answer==2) .* predChoice + (1-(answer==2)).* (1-predChoice);

    tempBeh = smoothdata(correct,'movmean',50); tempModel = smoothdata(predAcc,'movmean',50);
    rsq(i) = corr(tempBeh,tempModel)^2;

    %tempStim = stimulus; tempStim(stimulus==-1) = 2; tempResp = predChoice; tempResp(tempResp)
    %[tempModelBias, tempModel, ~, ~, ~] = fn_getAccBiasSmooth(stimulus, correct, 50);
    if strcmp(mouse{i},'zz063')
        fn_figureSmartDim('hSize',0.3,'widthHeightRatio',1); hold on;
        plot(tempBeh,'LineWidth',2)
        plot(tempModel,'LineWidth',2)
        disp(['R-square = ' num2str(corr(tempBeh,tempModel)^2,'%.3f')])
        xlim([0 2000]); xticks(0:500:2500); yticks(0.25:0.25:1);ylim([0.25 1])
    end
end


fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.4);
fn_plotComparison(rsq,'compType','boxplot'); ylim([0.5 1]); yticks(0.5:0.25:1)


%% -------------------VISUALIZE WEIGHT IN ALL MODELS------------------------
nPrev = 3;
mouse = {'zz054','zz062','zz063','zz066','zz067','zz068','zz069'};

% Action Bias visualization
wName = 'B';modelComparison = {'SB','SBA','SBRp','SBRn','SBRpRn'};
wIdx = {1,nPrev+1,nPrev+1,nPrev+1,nPrev*2+1};
w = loadW(mouse,modelComparison,wIdx,nPrev,true);
plotW_allModel(w,wName,modelComparison)

% Action history visualization - history 1
wName = 'A1';modelComparison = {'SA','SBA'}; wIdx = {1,1,1};
w = loadW(mouse,modelComparison,wIdx,nPrev,false);
plotW_allModel(w,wName,modelComparison)

% Action history visualization - history 2
wName = 'A2';modelComparison = {'SA','SBA'}; wIdx = {2,2,2};
w = loadW(mouse,modelComparison,wIdx,nPrev,false);
plotW_allModel(w,wName,modelComparison)

% Action history visualization - history 3
wName = 'A3';modelComparison = {'SA','SBA'}; wIdx = {3,3,3};
w = loadW(mouse,modelComparison,wIdx,nPrev,false);
plotW_allModel(w,wName,modelComparison)

% ActionReward history visualization - history 1
wName = 'Rp1';modelComparison = {'SRp','SRpRn','SBRp','SBRpRn'}; wIdx = {1,1,1,nPrev+1};
w = loadW(mouse,modelComparison,wIdx,nPrev,false);
plotW_allModel(w,wName,modelComparison)

% ActionReward history visualization - history 2
wName = 'Rp2';modelComparison = {'SRp','SRpRn','SBRp','SBRpRn'}; wIdx = {2,2,2,nPrev+2};
w = loadW(mouse,modelComparison,wIdx,nPrev,false);
plotW_allModel(w,wName,modelComparison)

% ActionReward history visualization - history 3
wName = 'Rp3';modelComparison = {'SRp','SRpRn','SBRp','SBRpRn'}; wIdx = {3,3,3,nPrev+3};
w = loadW(mouse,modelComparison,wIdx,nPrev,false);
plotW_allModel(w,wName,modelComparison)

%% -------------------VISUALIZE BIAS WEIGHT CONTRIBUTION TO CHOICE------------------------
wName = 'B';modelComparison = {'SBRpRn'}; wIdx = {nPrev*2+1};
wAccuracy = {};
figure; subplot(1,3,1); hold on;
w = loadW(mouse,modelComparison,wIdx,nPrev,false);w = w{1}; wChoice = fn_logistic(w);
plot((wChoice'),'Color',[0.8 0.8 0.8]); plot(nanmean((wChoice),1),'Color',[0 0 0])
xlim([1 3500]); ylim([0 1])
subplot(1,3,2); hold on;
w = loadW(mouse,modelComparison,wIdx,nPrev,true);w = w{1}; wChoice = fn_logistic(w);
plot((wChoice'),'Color',[0.8 0.8 0.8]); plot(nanmean((wChoice),1),'Color',[0 0 0])
xlim([1 3500]); ylim([0.5 1])

w = loadW(mouse,modelComparison,wIdx,nPrev,false);w = w{1}; wChoice = fn_logistic(w);
for i = 1:size(w,1)
    load(['C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\trialData\' filesep mouse{i} '_nPrev5.mat']);
    wAccuracy{i} = (y==1) .* (1-wChoice(i,1:length(y)))' + (y==2) .* (wChoice(i,1:length(y)))';
end
wAccuracy = fn_cell2matFillNan(wAccuracy); wAccuracySmooth = smoothdata(wAccuracy,1,'movmean',100);
subplot(1,3,3); hold on;
plot((wAccuracySmooth),'Color',[0.8 0.8 0.8]); plot(nanmean((wAccuracySmooth),2),'Color',[0 0 0])
xlim([1 3500]); ylim([0.4 0.9])
%% -------------------VISUALIZE HISTORY WEIGHT CONTRIBUTION TO CHOICE------------------------
wName = 'Rp';model = 'SRp'; wIdx = 1:3;
[choice,choicePred,accuracy] = loadMultipleW(mouse,model,wIdx,'actionXposRewardH',nPrev);
figure; subplot(1,2,1); hold on; binChoicePred = smoothdata(choicePred,1,'movmean',50);
plot(binChoicePred,'Color',[0.8 0.8 0.8]); plot(nanmean((binChoicePred),2),'Color',[0 0 0])
xlim([1 3500]); ylim([0.4 0.6]); title('Choice Prediction')
subplot(1,2,2); hold on; binAcc = smoothdata(accuracy,1,'movmean',50);
plot(binAcc,'Color',[0.8 0.8 0.8]); plot(nanmean(binAcc,2),'Color',[0 0 0])
xlim([1 3500]); ylim([0.4 0.6]); title('Accuracy')

%% -------------------VISUALIZE HISTORY WEIGHT CONTRIBUTION TO CHOICE------------------------
wName = 'Rn';model = 'SRn'; wIdx = 1:3;
[choice,choicePred,accuracy] = loadMultipleW(mouse,model,wIdx,'actionXnegRewardH',nPrev);
figure; subplot(1,2,1); hold on; binChoicePred = smoothdata(choicePred,1,'movmean',50);
plot(binChoicePred,'Color',[0.8 0.8 0.8]); plot(nanmean((binChoicePred),2),'Color',[0 0 0])
xlim([1 3500]); ylim([0.4 0.7]); title('Choice Prediction')
subplot(1,2,2); hold on; binAcc = smoothdata(accuracy,1,'movmean',50);
plot(binAcc,'Color',[0.8 0.8 0.8]); plot(nanmean(binAcc,2),'Color',[0 0 0])
xlim([1 3500]); ylim([0.4 0.7]); title('Accuracy')

%% -------------------VISUALIZE ACTION HISTORY WEIGHT CONTRIBUTION TO CHOICE------------------------
wName = 'A';model = 'SA'; wIdx = 1:3;
[choice,choicePred,accuracy] = loadMultipleW(mouse,model,wIdx,'actionH',nPrev);
figure; subplot(1,2,1); hold on; binChoicePred = smoothdata(choicePred,1,'movmean',50);
%temp = abs(binChoicePred-0.5)+0.5;
plot(binChoicePred,'Color',[0.8 0.8 0.8]); plot(nanmean((binChoicePred),2),'Color',[0 0 0])
xlim([1 3500]); ylim([0.4 0.73]); title('Choice Prediction')
subplot(1,2,2); hold on; binAcc = smoothdata(accuracy,1,'movmean',50);
plot(binAcc,'Color',[0.8 0.8 0.8]); plot(nanmean(binAcc,2),'Color',[0 0 0])
xlim([1 3500]); ylim([0.4 0.75]); title('Accuracy')

%% -------------------VISUALIZE ACTION HISTORY WEIGHT CONTRIBUTION TO CHOICE------------------------
[choice,choicePred,accuracy] = loadFullW(mouse,nPrev);
figure; subplot(1,2,1); hold on; binChoicePred = smoothdata(choicePred,1,'movmean',50);
%temp = abs(binChoicePred-0.5)+0.5;
plot(binChoicePred,'Color',[0.8 0.8 0.8]); plot(nanmean((binChoicePred),2),'Color',[0 0 0])
xlim([1 3500]); ylim([0.1 0.9]); title('Choice Prediction')
subplot(1,2,2); hold on; binAcc = smoothdata(accuracy,1,'movmean',50);
plot(binAcc,'Color',[0.8 0.8 0.8]); plot(nanmean(binAcc,2),'Color',[0 0 0])
xlim([1 3500]); ylim([0.4 0.8]); title('Accuracy')

%% -------------------All FUNCTIONS------------------------
function w = loadW(mouse,modelComparison,wIdx,nPrev,absFlag)
global loadPath;
for i = 1:length(modelComparison)
    tempW = {};
    for j = 1:length(mouse)
        if strcmp(modelComparison{i},'S') || strcmp(modelComparison{i},'SB')
            load([loadPath filesep mouse{j} 'psytrack_' modelComparison{i} '.mat'])
        else
            load([loadPath filesep mouse{j} 'psytrack_' modelComparison{i} '_nPrev' int2str(nPrev) '.mat'])
        end
        if absFlag; tempW{j} = abs(wMode(wIdx{i},:));
        else ; tempW{j} = wMode(wIdx{i},:);
        end
    end
    w{i} = fn_cell2matFillNan(tempW);
end
end


function [choice, choicePred,accuracy] = loadMultipleW(mouse,model,wIdx,varName,nPrev)

global loadPath;
for j = 1:length(mouse)
    if strcmp(model,'S') || strcmp(model,'SB')
        load([loadPath filesep mouse{j} 'psytrack_' model '.mat'])
    else
        load([loadPath filesep mouse{j} 'psytrack_' model '_nPrev' int2str(nPrev) '.mat'])
    end
    
    load(['C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\trialData' filesep mouse{j} '_nPrev5.mat']);
    behavVar = eval(varName);
    choicePred{j} = fn_logistic(sum(wMode(wIdx,:)' .* behavVar(:,wIdx),2));
    choice{j} = y;
    accuracy{j} = (y==1) .* (1-choicePred{j}) + (y==2) .* (choicePred{j});
    
end

choice = fn_cell2matFillNan(choice);
choicePred = fn_cell2matFillNan(choicePred);
accuracy = fn_cell2matFillNan(accuracy);

end


function [choice, choicePred,accuracy] = loadFullW(mouse,nPrev)
model = 'SBARpRn';
global loadPath;
for j = 1:length(mouse)
    if strcmp(model,'S') || strcmp(model,'SB')
        load([loadPath filesep mouse{j} 'psytrack_' model '.mat'])
    else
        load([loadPath filesep mouse{j} 'psytrack_' model '_nPrev' int2str(nPrev) '.mat'])
    end
    
    load(['C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\trialData' filesep mouse{j} '_nPrev5.mat']);
    
    regressor = cat(2,actionH(:,1:nPrev), actionXnegRewardH(:,1:nPrev), actionXposRewardH(:,1:nPrev), ones(size(actionH,1),1));
    
    choicePred{j} = fn_logistic(sum(regressor.* wMode(1:end-1,:)',2));
    choice{j} = y;
    accuracy{j} = (y==1) .* (1-choicePred{j}) + (y==2) .* (choicePred{j});
    
end

choice = fn_cell2matFillNan(choice);
choicePred = fn_cell2matFillNan(choicePred);
accuracy = fn_cell2matFillNan(accuracy);

end


function plotW_allModel(w,wName,modelComparison)
figure; [nRow,nColumn] = fn_sqrtInt(length(w));
for i = 1:length(w)
    subplot(nRow,nColumn,i); hold on;
    plot(w{i}','Color',[0.8 0.8 0.8]); plot(nanmean(w{i},1),'Color',[0 0 0])
    xlim([0 3500]); title([wName ' - ' modelComparison{i}]); xlabel('Trials'); ylabel('Weight')
    if i == 1; ylimm = ylim; else; ylim(ylimm); end
end
end

