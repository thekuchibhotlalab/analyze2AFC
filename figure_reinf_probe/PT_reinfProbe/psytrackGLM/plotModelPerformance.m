clear;
dataPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData';
sep = '\';
plotLim = [1 2400];
mice = {'zz054','zz062','zz063','zz066','zz067','zz068','zz069'};

%% 1.1 - MAKE ACCURACY LEARNING CURVE ACROSS MICE
learningCurveBin = 50;
accuracyCell = {};
for i = 1:length(mice)
    mouse = mice{i};
    behavior = load([dataPath sep mouse '.mat']);
    accuracyCell{i} = smoothdata(behavior.correct,'movmean', learningCurveBin);
end
f = figure;accuracyPlot = plotBlackGrayCurve(accuracyCell,'xlimm',plotLim,'ylimm',[0.3 1],...
    'color1',matlabColors(1,0.3),'color2',matlabColors(1,0.9));
%xlabel('trials'); ylabel('accuracy');

%% 1.2 - MAKE REGRESSION LEARNING CURVE ACROSS MICE
weightCell = {};
for i = 1:length(mice)
    mouse = mice{i};
    model = load([dataPath sep mouse 'psytrack_SBAArArs.mat']);
    maxAcc = max(accuracyCell{i});
    minAcc = min(accuracyCell{i});
    weightCell{i} = model.wMode(5,:);
    weightCell{i} = (1+exp(-weightCell{i})).^(-1);
    %weightCell{i} = (weightCell{i}-min(weightCell{i})) / (max(weightCell{i})...
    %    - min(weightCell{i})) * (maxAcc-minAcc)+minAcc;
end
%f = figure; 
weightPlot = plotBlackGrayCurve(weightCell,'xlimm',plotLim,'ylimm',[0.3 1],...
    'color1',matlabColors(2,0.3),'color2',matlabColors(2,0.9));
xlabel('trials'); ylabel('accuracy');

%% 2.1 - COMPARE BIAS ACROSS MICE
biasBin = 50; bias = {};bias_abs = {};
for i = 1:length(mice)
    mouse = mice{i};
    behavior = load([dataPath sep mouse '.mat']);
    
    bias{i} = nan(size(behavior.correct));
    for j = biasBin:length(behavior.correct)
        tempCorrect = behavior.correct((j-biasBin+1):j);
        tempStim = behavior.stimulus((j-biasBin+1):j);
        acc_s1 = sum(tempCorrect & tempStim==-1) / sum(tempStim==-1);
        acc_s2 = sum(tempCorrect & tempStim==1) / sum(tempStim==1);
        bias{i}(j) = (acc_s2 - acc_s1) / (acc_s1 + acc_s2);
        bias_abs{i}(j) = abs(bias{i}(j));
    end
    bias_abs{i}(1:biasBin-1) = nan;
end
f = figure; subplot(2,1,1)
biasPlot = plotBlackGrayCurve(bias,'baseline',0,'ylimm',[-1 1],'xlimm',plotLim);
subplot(2,1,2)
bias_absPlot = plotBlackGrayCurve(bias_abs,'baseline',0,'ylimm',[0 1],'xlimm',plotLim);

%% 2.2 - REGRESSION BIAS ACROSS MICE
biasWCell = {}; biasW_absCell = {};
for i = 1:length(mice)
    mouse = mice{i};
    model = load([dataPath sep mouse 'psytrack_SBAArArs.mat']);
    maxAcc = max(accuracyCell{i});
    minAcc = min(accuracyCell{i});
    biasWCell{i} = model.wMode(4,:);
    biasWCell{i} = (1+exp(-biasWCell{i})).^(-1)*2-1;
    biasW_absCell{i} = abs(biasWCell{i});
end
f = figure; subplot(2,1,1)
biasWPlot = plotBlackGrayCurve(biasWCell,'xlimm',plotLim,'ylimm',[-1 1],'baseline',0);
xlabel('trials'); ylabel('accuracy');
subplot(2,1,2)
biasW_absPlot = plotBlackGrayCurve(biasW_absCell,'xlimm',plotLim,'ylimm',[0 1],'baseline',0);
xlabel('trials'); ylabel('accuracy');

%% 2.3 - COMPARE BEHAVIORAL BIAS AND MODEL BIAS
figure;
bias_absPlot = plotBlackGrayCurve(bias_abs,'baseline',0,'ylimm',[0 1],'xlimm',plotLim,...
    'color1',matlabColors(1,0.3),'color2',matlabColors(1,0.9));
biasW_absPlot = plotBlackGrayCurve(biasW_absCell,'xlimm',plotLim,'ylimm',[0 1],'baseline',0,...
    'color1',matlabColors(2,0.3),'color2',matlabColors(2,0.9));



figure; subplot(2,1,1); hold on; 
plot([1 length(biasPlot)],[0 0],'LineWidth',3,'Color',[0.8 0.8 0.8])
plot(biasPlot,'Color',matlabColors(1,0.9),'LineWidth',2);
plot(biasWPlot,'Color',matlabColors(2,0.9),'LineWidth',2); 
xlim(plotLim); ylim([-0.5 0.5])
subplot(2,1,2); hold on
plot([1 length(bias_absPlot)],[0 0],'LineWidth',3,'Color',[0.8 0.8 0.8])
plot(bias_absPlot,'Color',matlabColors(1,0.9),'LineWidth',2);
plot(biasW_absPlot,'Color',matlabColors(2,0.9),'LineWidth',2); 
xlim(plotLim); ylim([0 1])

%% 3 - COMPARE MODEL PERFORMANCE

modelNames = {'S','SBAArArs','SBAAr','SBA','SB','SA','SAr','SArs'};
norm_logli = zeros(length(mice),length(modelNames));
for i = 1:length(mice)
    mouse = mice{i};
    for j = 1:length(modelNames)
        model = load([dataPath sep mouse 'psytrack_' modelNames{j} '.mat']);
        norm_logli(i,j) = model.xval_logli / length(model.xval_pL);
    end
end

norm_logli_diffS = norm_logli - repmat(norm_logli(:,1),[1 size(norm_logli,2)]);
norm_logli_diffS(:,1) = [];
figure;
bar(mean(norm_logli_diffS,1),'EdgeColor',[1 1 1],'FaceColor',matlabColors(2,0.9)); hold on;
for i = 1:size(norm_logli_diffS,2)
    scatter(i*ones(size(norm_logli_diffS,1),1) ,norm_logli_diffS(:,i),10,[0 0 0],'filled');
end
ylimm = ylim;ylim([0 ylimm(2)])
xticklabels({'+B+A+Ar+Ars','+B+A+Ar','+B+A','+B','+A','+Ar','+Ars'})
xtickangle(45); ylabel('log-likelihood increase')

figure;
bar(-mean(norm_logli(:,1:2),1),'EdgeColor',[1 1 1],'FaceColor',matlabColors(2,0.9)); hold on;
for i = 1:2
    scatter(i*ones(size(norm_logli,1),1),-norm_logli(:,i),10,[0 0 0],'filled');
end
xticklabels({'Stim Only',['+Bias '  '& Trial History']})
xtickangle(25); ylabel('Neg log-likeli (lower=better)');


%% 4.1 - REGRESSION ACTION HISTORY ACROSS MICE
wCell = {}; w_absCell = {};
for i = 1:length(mice)
    mouse = mice{i};
    model = load([dataPath sep mouse 'psytrack_SA.mat']);
    maxAcc = max(accuracyCell{i});
    minAcc = min(accuracyCell{i});
    wCell{i} = model.wMode(1,:);
    wCell{i} = (1+exp(-wCell{i})).^(-1)*2-1;
    w_absCell{i} = abs(wCell{i});
end
f = figure; subplot(2,1,1)
biasWPlot = plotBlackGrayCurve(wCell,'xlimm',plotLim,'ylimm',[-0.2 0.8],'baseline',0);
xlabel('trials'); ylabel('Action History');
subplot(2,1,2)
biasW_absPlot = plotBlackGrayCurve(w_absCell,'xlimm',plotLim,'ylimm',[0 0.5],'baseline',0);
xlabel('trials'); ylabel('Abs Action History');

%% 4.2 - REGRESSION ACTION HISTORY ACROSS MICE
mice = {'zz062','zz067','zz068'};
wCell = {}; w_absCell = {};
for i = 1:length(mice)
    mouse = mice{i};
    model = load([dataPath sep mouse 'psytrack_SAr.mat']);
    maxAcc = max(accuracyCell{i});
    minAcc = min(accuracyCell{i});
    wCell{i} = model.wMode(1,:);
    wCell{i} = (1+exp(-wCell{i})).^(-1)*2-1;
    w_absCell{i} = abs(wCell{i});
end
f = figure; subplot(2,1,1)
biasWPlot = plotBlackGrayCurve(wCell,'xlimm',plotLim,'ylimm',[-0.2 0.8],'baseline',0);
xlabel('trials'); ylabel('Action X Reward');
subplot(2,1,2)
biasW_absPlot = plotBlackGrayCurve(w_absCell,'xlimm',plotLim,'ylimm',[0 0.5],'baseline',0);
xlabel('trials'); ylabel('Abs Action X Reward');

%% 4.3 - REGRESSION ACTION HISTORY ACROSS MICE
mice = {'zz062','zz067','zz068'};
wCell = {}; w_absCell = {};
for i = 1:length(mice)
    mouse = mice{i};
    model = load([dataPath sep mouse 'psytrack_SArs.mat']);
    maxAcc = max(accuracyCell{i});
    minAcc = min(accuracyCell{i});
    wCell{i} = model.wMode(1,:);
    wCell{i} = (1+exp(-wCell{i})).^(-1)*2-1;
    w_absCell{i} = abs(wCell{i});
end
f = figure; subplot(2,1,1)
biasWPlot = plotBlackGrayCurve(wCell,'xlimm',plotLim,'ylimm',[-0.2 0.8],'baseline',0);
xlabel('trials'); ylabel('Action X Reward X Stim');
subplot(2,1,2)
biasW_absPlot = plotBlackGrayCurve(w_absCell,'xlimm',plotLim,'ylimm',[0 0.7],'baseline',0);
xlabel('trials'); ylabel('Abs Action X Reward X Stim');


%% 5.1 - ACTION HISTORY EFFECT IN BEHAVIOR
mice = {'zz054','zz062','zz063','zz066','zz067','zz068','zz069'};
binSize = 300; 
binWindow = 10;
learningBin = {}; historyEffectVec = {}; historyEffect = {};
for i = 1:length(mice)
    mouse = mice{i};
    behavior = load([dataPath sep mouse '.mat']);
    
    learningBinStart = binSize:binWindow:length(behavior.stimulus);
    learningBin{i} = learningBinStart - binSize/2;
    
    for j = 1:length(learningBinStart)
        tempIdx = learningBinStart(j)-binSize+1 : learningBinStart(j);
        [val,vec] = getTrialHistory(behavior,tempIdx,behavior.actionH);
        historyEffectVec{i}(:,j) = vec;
        historyEffect{i}(j) = val; 
    end
end

f = figure;accuracyPlot = plotBlackGrayCurve(historyEffect,'xlimm',[0 250],'baseline',0,'ylim',[-0.2 0.4]);
xlabel('trials'); ylabel('delta p(choice)');
xticks(50:50:250);xticklabels({'500','1000','1500','2000','2500'})


%% 5.2 - ACTION HISTORY EFFECT IN BEHAVIOR
binSize = 300; 
binWindow = 10;
learningBin = {}; historyEffectVec = {}; historyEffect = {};
for i = 1:length(mice)
    mouse = mice{i};
    behavior = load([dataPath sep mouse '.mat']);
    
    learningBinStart = binSize:binWindow:length(behavior.stimulus);
    learningBin{i} = learningBinStart - binSize/2;
    
    for j = 1:length(learningBinStart)
        tempIdx = learningBinStart(j)-binSize+1 : learningBinStart(j);
        [val,vec] = getTrialHistory(behavior,tempIdx,behavior.actionXrewardH);
        historyEffectVec{i}(:,j) = vec;
        historyEffect{i}(j) = val; 
    end
end

f = figure;accuracyPlot = plotBlackGrayCurve(historyEffect,'xlimm',[0 250],'baseline',0,'ylim',[-0.2 0.4]);
xlabel('trials'); ylabel('delta p(choice)');
xticks(50:50:250);xticklabels({'500','1000','1500','2000','2500'})

%% ALL FUNCTIONS
function avgPlot = plotBlackGrayCurve(matCell,varargin)

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('baseline', 0.5)
p.addParameter('ylimm', [0 1])
p.addParameter('xlimm', [])
p.addParameter('color1', [0.8 0.8 0.8])
p.addParameter('color2', [0.1 0.1 0.1])

p.parse(varargin{:});

baseline = p.Results.baseline; 
ylimm = p.Results.ylimm;
xlimm = p.Results.xlimm;
color1 = p.Results.color1;
color2 = p.Results.color2;

tempMat = nan(length(matCell),10000);
for i = 1:length(matCell)
    plot(matCell{i},'Color',color1);hold on;
    tempMat(i,1:length(matCell{i})) = matCell{i};
end
avgPlot = nanmean(tempMat,1); avgPlot(isnan(avgPlot)) = [];
plot([1 length(avgPlot)],[baseline baseline],'LineWidth',3,'Color',[0.8 0.8 0.8]) 
plot(avgPlot,'Color',color2,'LineWidth',2);
if isempty(xlimm); xlim([1 length(avgPlot)]); else; xlim(xlimm); end
ylim(ylimm)

end

function [val,vec]=getTrialHistory(behavior,tempIdx,history)

ah1Flag = history(tempIdx,1) == -1; ah2Flag = history(tempIdx,1) == 1; 
s1Flag = behavior.stimulus(tempIdx)==-1; s2Flag = behavior.stimulus(tempIdx)==1;
action = behavior.y(tempIdx);

% P(action1) given AH == 1
% Balance current S1 and current S2 proportion
tempS1 = action(ah1Flag & s1Flag); tempS2 = action(ah1Flag & s2Flag); 
ah1Action1 = [sum(tempS1==1)/length(tempS1)  sum(tempS2==1)/length(tempS2)];
% P(action1) given AH == 2
tempS1 = action(ah2Flag & s1Flag); tempS2 = action(ah2Flag & s2Flag); 
ah2Action1 = [sum(tempS1==1)/length(tempS1)  sum(tempS2==1)/length(tempS2)];

% P(action2) given AH == 1
tempS1 = action(ah1Flag & s1Flag); tempS2 = action(ah1Flag & s2Flag); 
ah1Action2 = [sum(tempS1==2)/length(tempS1)  sum(tempS2==2)/length(tempS2)];
% P(action2) given AH == 2
tempS1 = action(ah2Flag & s1Flag); tempS2 = action(ah2Flag & s2Flag); 
ah2Action2 = [sum(tempS1==2)/length(tempS1)  sum(tempS2==2)/length(tempS2)];

vec = [mean(ah1Action1) mean(ah2Action1) mean(ah1Action2) mean(ah2Action2)];
val = (mean(ah1Action1)-mean(ah2Action1)-mean(ah1Action2)+mean(ah2Action2))/2;


end

