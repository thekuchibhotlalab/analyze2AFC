%% Plot example animal
clear; global loadPath;
loadPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\psyTrackFit\';
mouse ='zzRandom'; modelComparison = 'SB'; iter = 17;
load([loadPath filesep mouse 'psytrack_' modelComparison int2str(iter) '.mat']);
sW = (wMode(2,:)); bW = (wMode(1,:)); choice(choice==2) =-1;
acc = smoothdata(stimulus==choice','movmean',100); bias = smoothdata(choice,'movmean',100)/2;
swPred = fn_logistic(sW); bwPred = 2*(fn_logistic(bW)-0.5);

fn_figureSmartDim('hSize',0.3,'widthHeightRatio',0.9); hold on;
plot(acc,'Color', matlabColors(1,0.9), 'LineWidth',2)
plot(swPred,'Color', matlabColors(2,0.9), 'LineWidth',2);
ylim([0 1]); yticks(0:0.5:1); 
xlim([0 2000]);xticks(0:500:2000)
fn_figureSmartDim('hSize',0.3,'widthHeightRatio',0.9); hold on;
plot(bias,'Color', matlabColors(1,0.9), 'LineWidth',2)
plot(bwPred,'Color', matlabColors(2,0.9), 'LineWidth',2);
ylim([-0.5 0.5]); yticks(-0.5:0.5:0.5); 
xlim([0 2000]);xticks(0:500:2000)


%% Plot average weights
clear;
global loadPath;
loadPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\psyTrackFit\';

mouse ='zzRandom';
modelComparison = 'SB';
learningLim = [];
iter = 0:18;
sW = {}; bW = {};
for i = 1:length(iter)
    load([loadPath filesep mouse 'psytrack_' modelComparison int2str(iter(i)) '.mat']);
    sW{i} = (wMode(2,:)); bW{i} = (wMode(1,:)); 
end
sWRaw = fn_cell2matFillNan(sW); bWRaw = fn_cell2matFillNan(bW);
sW = abs(fn_cell2matFillNan(sW)); bW = abs(fn_cell2matFillNan(bW));

fn_figureSmartDim('hSize',0.3,'widthHeightRatio',0.9); hold on;
f_errorbar = fn_plotFillErrorbar(1:size(sW,2),nanmean(sW,1),nanstd(sW,0,1)./sqrt(sum(~isnan(sW),1)),...
    matlabColors(1),'faceAlpha',0.2,'LineStyle','none');
plot(nanmean(sW,1),'Color', matlabColors(1,0.9), 'LineWidth',2);

f_errorbar = fn_plotFillErrorbar(1:size(bW,2),nanmean(bW,1),nanstd(bW,0,1)./sqrt(sum(~isnan(bW),1)),...
    matlabColors(2),'faceAlpha',0.2,'LineStyle','none');
plot(nanmean(bW,1),'Color', matlabColors(2,0.9), 'LineWidth',2);
ylim([0 2])
% PLOT the average weights
sWMean = nanmean(sW,2);
bWMean = nanmean(bW,2);

tempPlot = {sWMean,bWMean};
fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.5);hold on;
plot([0 6],[0 0],'LineWidth',2,'Color',[0.8 0.8 0.8])
violin(tempPlot,'facecolor',matlabColors([1:4 7]),'edgecolor',[0.4 0.4 0.4],'mc',[0.4 0.4 0.4],'medc',[],'facealpha',0.5);

for i = 1:length(tempPlot)
    scatter(ones(1,length(tempPlot{i}))*i,tempPlot{i},8,[0.4 0.4 0.4],'filled');
end
xlim([0.5 2.5]); xticks(1:2); xticklabels({'s';'b'}); ylim([-0.5 2.5])
