function plotLearningCurveTrialAll(mice,selectProtocol)

if nargin == 1; selectProtocol = {'puretone'}; end
global sep;

rootPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\';
loadPath = [rootPath sep 'trialData\'];
figPath = [rootPath sep 'figure\' sep 'learningCurveTrialAll' sep]; mkdir(figPath);
trialData = cellfun(@(x)fn_selectProtocol(loadPath, x, selectProtocol),mice,'UniformOutput',false);


%[~, accuracy,~,~,~,~,~,~,~,~,~,~,probeData,~, ~,~,probeTrialNum,~,~]
outmat = cellfun(@(x)alignDaysKeepMiss(x),trialData,'UniformOutput',false);

accuracy = cellfun(@(x)(x.accuracy),outmat,'UniformOutput',false);
probeData = cellfun(@(x)(x.probeData),outmat,'UniformOutput',false);
probeTrialNum = cellfun(@(x)(x.probeTrialNum),outmat,'UniformOutput',false);

trialLim = 2400;
trialBin = 0:300:trialLim;
probeAllAni = nan(length(probeData),length(trialBin)-1);

for i= 1:length(probeData)
    trialNum = probeTrialNum{i}; 
    trialLimFlag = trialNum>trialLim; trialNum(trialLimFlag) = [];
    tempPerf = probeData{i}(~trialLimFlag,:);
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
accuracy = cellfun(@(x)(x(1:trialLim)),accuracy,'UniformOutput',false);
accuracy = fn_cell2mat(accuracy,2);

% PLOT 1
figure;
hold on;
f_errorbar = fn_plotFillErrorbar(trialBinPlot(~nanflag),meanProbePerf(~nanflag),semProbePerf(~nanflag),...
    matlabColors(2),'faceAlpha',0.2,'LineStyle','none');
f_errorbar = fn_plotFillErrorbar(1:size(accuracy,1),(mean(accuracy,2))',...
    (std(accuracy,0,2)./sqrt(size(accuracy,2)))',...
    matlabColors(1),'faceAlpha',0.2,'LineStyle','none');
% PLOT REINF
plot(mean(accuracy,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
% PLOT PROBE AGAIN
plot(trialBinPlot(~nanflag),meanProbePerf(~nanflag),'--o','Color',matlabColors(2,0.9),'LineWidth',2)

% PLOT 2
figure;
hold on;
% PLOT REINF
for j = 1:size(accuracy,2)
    plot(accuracy(:,j),'Color', matlabColors(1,0.2), 'LineWidth',1);
end
for j = 1:size(probeAllAni,1)
    plot(trialBinPlot,probeAllAni(j,:),'--o','Color', matlabColors(2,0.4), 'LineWidth',1)
end
plot(mean(accuracy,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
% PLOT PROBE AGAIN
plot(trialBinPlot(~nanflag),meanProbePerf(~nanflag),'--o','Color',matlabColors(2,0.9),'LineWidth',2)

end