%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);

probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
%% PLOT 1.0 -- Subsample within DAY
probeDay = mouseMega.objFun('computeProbeAlignByTrial',{'day',[reinfThre probeThre],probeTrialBin});
probeDayPerf = mouseMega.objFun('computeProbeAlignByPerf',{'day',[0.6 0.8],400});
% cell identity -- [probeDataAllTrial,probeByDay,nProbePerDay,probeByDayBin,trialLim]
%% PLOT 1.1 -- reinf vs. probe comparison
plotSubsampleDay(probeDay,'Select by aligned trial, ');
%% PLOT 1.1 -- reinf vs. probe comparison, by performance
plotSubsampleDay(probeDayPerf,'Select by performance, ');
%% PLOT 2.0 -- Subsample within Bef & Aft
probeBefAft = mouseMega.objFun('computeProbeAlignByTrial',{'befAft',[reinfThre probeThre],probeTrialBin});
probeBefAftPerf = mouseMega.objFun('computeProbeAlignByPerf',{'befAft',[0.6 0.8],400});
%% PLOT 2.1 -- reinf vs. probe comparison
plotBefAftCombined(probeBefAft, 'Select by aligned trial, ');
%% PLOT 2.2 -- bef vs. aft comparison
plotBefAft(probeBefAft);
%% PLOT 2.3 -- bef vs. aft comparison, by performance
plotBefAft(probeBefAftPerf);
%%

function plotBefAft(probeCount, addTitle)
if nargin == 1; addTitle = ''; end
% plot by day
probeByDay = fn_cell2mat(cellfun(@(x)((abs(x.probe))),probeCount{2},'UniformOutput',false),1); 
befByDay = fn_cell2mat(cellfun(@(x)((abs(x.bef))),probeCount{2},'UniformOutput',false),1); 
aftByDay = fn_cell2mat(cellfun(@(x)((abs(x.aft))),probeCount{2},'UniformOutput',false),1); 

figure;
subplot(2,1,1);
fn_boxplot({befByDay(:,1),probeByDay(:,1),aftByDay(:,1)},'paired',true); xticklabels({'bef','probe','aft'})

subplot(2,1,2);
fn_boxplot({befByDay(:,2),probeByDay(:,2),aftByDay(:,2)},'paired',true); xticklabels({'bef','probe','aft'})
[p,tbl,stats] = anova1([(befByDay(:,1)) (probeByDay(:,1)) (aftByDay(:,1))]);
d = multcompare(stats);
[p,tbl,stats] = anova1([(befByDay(:,2)) (probeByDay(:,2)) (aftByDay(:,2))]);
d = multcompare(stats);
% plot by animal
probeAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,1)'))),probeCount{2},'UniformOutput',false)); 
probeBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,2)'))),probeCount{2},'UniformOutput',false)); 
befAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.bef(:,1)'))),probeCount{2},'UniformOutput',false)); 
befBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.bef(:,2)'))),probeCount{2},'UniformOutput',false)); 
aftAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.aft(:,1)'))),probeCount{2},'UniformOutput',false)); 
aftBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.aft(:,2)'))),probeCount{2},'UniformOutput',false)); 

figure;
subplot(2,1,1);
fn_boxplot({nanmean(befAcc,2),nanmean(probeAcc,2),nanmean(aftAcc,2)},'paired',true); xticklabels({'bef','probe','aft'})
subplot(2,1,2);
fn_boxplot({nanmean(befBias,2),nanmean(probeBias,2),nanmean(aftBias,2)},'paired',true); xticklabels({'bef','probe','aft'})
[p,tbl,stats] = anova1([nanmean(befAcc,2) nanmean(probeAcc,2) nanmean(aftAcc,2)]);
d = multcompare(stats);
[p,tbl,stats] = anova1([nanmean(befBias,2),nanmean(probeBias,2),nanmean(aftBias,2)]);
d = multcompare(stats);


end

function plotBefAftCombined(probeCount, addTitle)
if nargin == 1; addTitle = ''; end
% plot by day
probeByDay = fn_cell2mat(cellfun(@(x)((abs(x.probe))),probeCount{2},'UniformOutput',false),1); 
befByDay = fn_cell2mat(cellfun(@(x)((abs(x.bef))),probeCount{2},'UniformOutput',false),1); 
aftByDay = fn_cell2mat(cellfun(@(x)((abs(x.aft))),probeCount{2},'UniformOutput',false),1); 
reinfByDay = nanmean(cat(3,befByDay,aftByDay),3);
t = fn_barPlotReinfProbe(reinfByDay(:,1),reinfByDay(:,2),probeByDay(:,1),probeByDay(:,2),'ttest'); 
subplot(2,1,1);title([addTitle 'Comp by day, ttest' newline t{1}])
t = fn_barPlotReinfProbe(reinfByDay(:,1),reinfByDay(:,2),probeByDay(:,1),probeByDay(:,2),'signrank'); 
subplot(2,1,1);title([addTitle 'Comp by day, signrank' newline t{1}])

% plot by animal
probeAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,1)'))),probeCount{2},'UniformOutput',false)); 
probeBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,2)'))),probeCount{2},'UniformOutput',false)); 
befAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.bef(:,1)'))),probeCount{2},'UniformOutput',false)); 
befBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.bef(:,2)'))),probeCount{2},'UniformOutput',false)); 
aftAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.aft(:,1)'))),probeCount{2},'UniformOutput',false)); 
aftBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.aft(:,2)'))),probeCount{2},'UniformOutput',false)); 

reinfAcc = nanmean(cat(3,befAcc,aftAcc),3); reinfBias = nanmean(cat(3,befBias,aftBias),3);

t = fn_barPlotReinfProbe(nanmean(reinfAcc,2),nanmean(reinfBias,2),nanmean(probeAcc,2),nanmean(probeBias,2),'ttest');
subplot(2,1,1);title([addTitle 'Comp by animal, ttest' newline t{1}])
t = fn_barPlotReinfProbe(nanmean(reinfAcc,2),nanmean(reinfBias,2),nanmean(probeAcc,2),nanmean(probeBias,2),'signrank');
subplot(2,1,1);title([addTitle 'Comp by animal, signrank' newline t{1}])
end

function plotSubsampleDay(probeCount, addTitle)
if nargin == 1; addTitle = ''; end
% plot by day
probeByDay = fn_cell2mat(cellfun(@(x)((abs(x.probe))),probeCount{2},'UniformOutput',false),1); 
reinfByDay = fn_cell2mat(cellfun(@(x)((abs(x.reinf))),probeCount{2},'UniformOutput',false),1); 
t = fn_barPlotReinfProbe(reinfByDay(:,1),reinfByDay(:,2),probeByDay(:,1),probeByDay(:,2),'ttest'); 
subplot(2,1,1);title([addTitle 'Comp by day, ttest' newline t{1}])
t = fn_barPlotReinfProbe(reinfByDay(:,1),reinfByDay(:,2),probeByDay(:,1),probeByDay(:,2),'signrank'); 
subplot(2,1,1);title([addTitle 'Comp by day, signrank' newline t{1}])

% plot by animal
probeAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,1)'))),probeCount{2},'UniformOutput',false)); 
probeBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,2)'))),probeCount{2},'UniformOutput',false)); 
reinfAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.reinf(:,1)'))),probeCount{2},'UniformOutput',false)); 
reinfBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.reinf(:,2)'))),probeCount{2},'UniformOutput',false)); 

t = fn_barPlotReinfProbe(nanmean(reinfAcc,2),nanmean(reinfBias,2),nanmean(probeAcc,2),nanmean(probeBias,2),'ttest');
subplot(2,1,1);title([addTitle 'Comp by animal, ttest' newline t{1}])
t = fn_barPlotReinfProbe(nanmean(reinfAcc,2),nanmean(reinfBias,2),nanmean(probeAcc,2),nanmean(probeBias,2),'signrank');
subplot(2,1,1);title([addTitle 'Comp by animal, signrank' newline t{1}])
end

function t = fn_barPlotReinfProbe(reinfAcc,reinfBias,probeAcc,probeBias,test)
figure;
subplot(2,1,1);fn_boxplot({reinfAcc,probeAcc},'paired',true,'dotLoc','side',...
    'boxplotArgIn',{}); 
xticks([1 2]); xticklabels({'reinf','probe'}); ylim([0.25 1.0]); yticks([0.25 0.5 0.75 1.0])
if strcmp(test,'ttest')
    [h,p] = ttest(reinfAcc,probeAcc);  t{1} = ['accuracy, p= ' num2str(p,'%.4f')] ;
else
    [p,h] = signrank(reinfAcc,probeAcc); t{1} = ['accuracy, p= ' num2str(p,'%.4f')];
end
title(t{1});

subplot(2,1,2);fn_boxplot({reinfBias,probeBias},'paired',true,'dotLoc','side',...
    'boxplotArgIn',{}); 
xticks([1 2]); xticklabels({'reinf','probe'}); ylim([0 1.0]); yticks([0 0.25 0.5 0.75 1.0])
if strcmp(test,'ttest')
    [h,p] = ttest(reinfBias,probeBias);  t{2} = ['choice bias, p= ' num2str(p,'%.4f')];
else
    [p,h] = signrank(reinfBias,probeBias);  t{2} = ['choice bias, p= ' num2str(p,'%.4f')];
end
title(t{2});


end


