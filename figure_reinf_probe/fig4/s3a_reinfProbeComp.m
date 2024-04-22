%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);

probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;

%% PLOT a -- Early learning, all sessions
probeDay = mouseMega.objFun('computeProbeAlignByTrial',{'day',[reinfThre probeThre],probeTrialBin});
%%
probeAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,1)'))),probeDay{2},'UniformOutput',false)); 
probeBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,2)'))),probeDay{2},'UniformOutput',false)); 
reinfAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.reinf(:,1)'))),probeDay{2},'UniformOutput',false)); 
reinfBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.reinf(:,2)'))),probeDay{2},'UniformOutput',false)); 
figure;
plot([0 1],[0 0],'LineWidth',0.5,'Color',[0.8 0.8 0.8]);
[r2,p2] = corrcoef(reinfBias(:),probeBias(:),'rows','complete');
scatter(reinfBias(:),probeBias(:),30,matlabColors(1,0.9),'filled'); hold on;
lm = fitlm(reinfBias(:),probeBias(:));
tempX = [0 0.85]; tempY = tempX * lm.Coefficients{2,1} + lm.Coefficients{1,1};

plot(tempX,tempY,'Color',matlabColors(1,0.9),'LineWidth',0.5);
tempAxis = 0:0.01:0.85; [Ypred,YCI] = predict(lm, tempAxis');
fill([tempAxis fliplr(tempAxis)], [YCI(:,1);flipud(YCI(:,2))]',matlabColors(1,0.9),'LineStyle','None','faceAlpha',0.2);
title(['corr = ' num2str(r2(1,2)) ', p = ' num2str(p2(1,2)) ', slope = ' num2str(lm.Coefficients.Estimate(2))])
%% PLOT b-- All learning, subsample within DAY,
probeDay = mouseMega.objFun('computeProbeAlignByTrial',{'befAft',[reinfThre probeThre],probeTrialBin});
%% PLOT 1.1 -- reinf vs. probe comparison
plotSubsampleAnimal(probeDay);

%% PLOT 1.2 -- reinf vs. probe comparison
plotBefAft(probeDay);

%% PLOT 1.2 -- reinf vs. probe bias scatter

befAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.bef(:,1)'))),probeDay{2},'UniformOutput',false)); 
befBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.bef(:,2)'))),probeDay{2},'UniformOutput',false)); 
aftAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.aft(:,1)'))),probeDay{2},'UniformOutput',false)); 
aftBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.aft(:,2)'))),probeDay{2},'UniformOutput',false)); 
probeAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,1)'))),probeDay{2},'UniformOutput',false)); 
probeBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,2)'))),probeDay{2},'UniformOutput',false)); 
figure; hold on; 
plot([0 1],[0 1],'LineWidth',0.5,'Color',[0.8 0.8 0.8]);

[r1,p1] = corrcoef(befBias(:),probeBias(:),'rows','complete');
[r2,p2] = corrcoef(aftBias(:),probeBias(:),'rows','complete');
scatter(befBias(:),probeBias(:),30,matlabColors(1),'filled'); 
scatter(aftBias(:),probeBias(:),30,matlabColors(3),'filled'); 
lm = fitlm(befBias(:),probeBias(:)); lm2 = fitlm(aftBias(:),probeBias(:));
tempX = [0 1]; tempY = tempX * lm.Coefficients{2,1} + lm.Coefficients{1,1};
plot(tempX,tempY,'Color',matlabColors(1,0.9),'LineWidth',0.5);
tempX = [0 1]; tempY = tempX * lm2.Coefficients{2,1} + lm2.Coefficients{1,1};
plot(tempX,tempY,'Color',matlabColors(3,0.9),'LineWidth',0.5);

tempAxis = 0:0.01:1; [Ypred,YCI] = predict(lm, tempAxis');
fill([tempAxis fliplr(tempAxis)], [YCI(:,1);flipud(YCI(:,2))]',matlabColors(1,0.9),'LineStyle','None','faceAlpha',0.2);
[Ypred2,YCI2] = predict(lm2, tempAxis');
fill([tempAxis fliplr(tempAxis)], [YCI2(:,1);flipud(YCI2(:,2))]',matlabColors(3,0.9),'LineStyle','None','faceAlpha',0.2);
title(['corr = ' num2str(r1(1,2)) ', p = ' num2str(p1(1,2)) ', slope = ' num2str(lm.Coefficients.Estimate(2)) newline...
    'corr = ' num2str(r2(1,2)) ', p = ' num2str(p2(1,2)) ', slope = ' num2str(lm2.Coefficients.Estimate(2))]) 
%% PLOT 1.2 -- reinf vs. probe bias scatter
befAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.bef(:,1)'))),probeDay{2},'UniformOutput',false)); 
befBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.bef(:,2)'))),probeDay{2},'UniformOutput',false)); 
aftAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.aft(:,1)'))),probeDay{2},'UniformOutput',false)); 
aftBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.aft(:,2)'))),probeDay{2},'UniformOutput',false)); 
probeAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,1)'))),probeDay{2},'UniformOutput',false)); 
probeBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,2)'))),probeDay{2},'UniformOutput',false)); 

figure; hold on; 
[r1,p1] = corrcoef(befBias(:),befBias(:)-probeBias(:),'rows','complete');
lm = fitlm(befBias(:),befBias(:)-probeBias(:));
scatter(befBias(:),befBias(:)-probeBias(:),30,matlabColors(1),'filled'); 
title(['corr = ' num2str(r1(1,2)) ', p = ' num2str(p1(1,2)) ', slope = ' num2str(lm.Coefficients.Estimate(2))])

figure; hold on; 
[r1,p1] = corrcoef(befBias(:),probeAcc(:)-befAcc(:),'rows','complete');
[r2,p2] = corrcoef(aftBias(:),probeAcc(:)-aftAcc(:),'rows','complete');
scatter(befBias(:),probeAcc(:)-befAcc(:),30,matlabColors(1),'filled'); 
scatter(aftBias(:),probeAcc(:)-aftAcc(:),30,matlabColors(3),'filled'); 
lm = fitlm(befBias(:),probeAcc(:)-befAcc(:)); lm2 = fitlm(aftBias(:),probeAcc(:)-aftAcc(:));
tempX = [0 1]; tempY = tempX * lm.Coefficients{2,1} + lm.Coefficients{1,1};
plot(tempX,tempY,'Color',matlabColors(1,0.9),'LineWidth',0.5);
tempX = [0 1]; tempY = tempX * lm2.Coefficients{2,1} + lm2.Coefficients{1,1};
plot(tempX,tempY,'Color',matlabColors(3,0.9),'LineWidth',0.5);

tempAxis = 0:0.01:1; [Ypred,YCI] = predict(lm, tempAxis');
fill([tempAxis fliplr(tempAxis)], [YCI(:,1);flipud(YCI(:,2))]',matlabColors(1,0.9),'LineStyle','None','faceAlpha',0.2);
[Ypred2,YCI2] = predict(lm2, tempAxis');
fill([tempAxis fliplr(tempAxis)], [YCI2(:,1);flipud(YCI2(:,2))]',matlabColors(3,0.9),'LineStyle','None','faceAlpha',0.2);
title(['corr = ' num2str(r1(1,2)) ', p = ' num2str(p1(1,2)) ', slope = ' num2str(lm.Coefficients.Estimate(2)) newline...
    'corr = ' num2str(r2(1,2)) ', p = ' num2str(p2(1,2)) ', slope = ' num2str(lm2.Coefficients.Estimate(2))]) 
ylim([-0.5 0.5]); xlim([0 1])

% plot correlation of probe bias and reduction of bias
figure; hold on; 
[r1,p1] = corrcoef(befBias(:),probeBias(:)-befBias(:),'rows','complete');
[r2,p2] = corrcoef(aftBias(:),probeBias(:)-aftBias(:),'rows','complete');
scatter(befBias(:),probeBias(:)-befBias(:),30,matlabColors(1),'filled'); 
scatter(aftBias(:),probeBias(:)-aftBias(:),30,matlabColors(3),'filled'); 
lm = fitlm(befBias(:),probeBias(:)-befBias(:)); lm2 = fitlm(aftBias(:),probeBias(:)-aftBias(:));
tempX = [0 1]; tempY = tempX * lm.Coefficients{2,1} + lm.Coefficients{1,1};
plot(tempX,tempY,'Color',matlabColors(1,0.9),'LineWidth',0.5);
tempX = [0 1]; tempY = tempX * lm2.Coefficients{2,1} + lm2.Coefficients{1,1};
plot(tempX,tempY,'Color',matlabColors(3,0.9),'LineWidth',0.5);

tempAxis = 0:0.01:1; [Ypred,YCI] = predict(lm, tempAxis');
fill([tempAxis fliplr(tempAxis)], [YCI(:,1);flipud(YCI(:,2))]',matlabColors(1,0.9),'LineStyle','None','faceAlpha',0.2);
[Ypred2,YCI2] = predict(lm2, tempAxis');
fill([tempAxis fliplr(tempAxis)], [YCI2(:,1);flipud(YCI2(:,2))]',matlabColors(3,0.9),'LineStyle','None','faceAlpha',0.2);
title(['corr = ' num2str(r1(1,2)) ', p = ' num2str(p1(1,2)) ', slope = ' num2str(lm.Coefficients.Estimate(2)) newline...
    'corr = ' num2str(r2(1,2)) ', p = ' num2str(p2(1,2)) ', slope = ' num2str(lm2.Coefficients.Estimate(2))]) 
ylim([-0.8 0.8]); xlim([0 1])
%% FUNCTIONS
function plotSubsampleAnimal(probeCount)
% plot by animal
probeAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,1)'))),probeCount{2},'UniformOutput',false)); 
probeBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,2)'))),probeCount{2},'UniformOutput',false)); 
befAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.bef(:,1)'))),probeCount{2},'UniformOutput',false)); 
befBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.bef(:,2)'))),probeCount{2},'UniformOutput',false)); 
aftAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.aft(:,1)'))),probeCount{2},'UniformOutput',false)); 
aftBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.aft(:,2)'))),probeCount{2},'UniformOutput',false)); 

reinfAcc = nanmean(cat(3,befAcc,aftAcc),3); reinfBias = nanmean(cat(3,befBias,aftBias),3);

fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.4); 
subplot(2,1,1);hold on;
p = fn_plotComparison({nanmean(reinfAcc,2),nanmean(probeAcc,2)},'compType', 'errorbarWithDot','paired',true,'test','signrank'); 
title(['Accuracy, p=' num2str(p)]); ylim([0.5 1]);  xlim([0.6 2.4]); yticks([0.5 0.75 1.0]);

subplot(2,1,2); hold on;
p = fn_plotComparison({nanmean(reinfBias,2),nanmean(probeBias,2)},'compType', 'errorbarWithDot','paired',true,'test','signrank'); 
title(['Choice bias, p=' num2str(p)]); ylim([0 0.6]); yticks(0:0.2:0.6); xlim([0.6 2.4])


% plot by day
probeByDay = fn_cell2mat(cellfun(@(x)((abs(x.probe))),probeCount{2},'UniformOutput',false),1); 
befByDay = fn_cell2mat(cellfun(@(x)((abs(x.bef))),probeCount{2},'UniformOutput',false),1); 
aftByDay = fn_cell2mat(cellfun(@(x)((abs(x.aft))),probeCount{2},'UniformOutput',false),1); 
reinfByDay = nanmean(cat(3,befByDay,aftByDay),3);

% bar plot
fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.4); 
subplot(2,1,1);
p = fn_plotComparison({reinfByDay(:,1),probeByDay(:,1)},'compType', 'errorbarWithDot','paired',true,'test','signrank');
title(['Accuracy, p=' num2str(p)]); ylim([0.25 1]);  xlim([0.6 2.4]); yticks([0.25 0.5 0.75 1.0]);

subplot(2,1,2);
p = fn_plotComparison({reinfByDay(:,2),probeByDay(:,2)},'compType', 'errorbarWithDot','paired',true,'test','signrank'); 
title(['Choice bias, p=' num2str(p)]); ylim([0 1.0]); yticks(0:0.5:1.0); xlim([0.6 2.4])

% scatter plot
fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.5); 
subplot(2,1,1); hold on;
scatter(reinfByDay(:,1),probeByDay(:,1),10,'filled'); plot([0 1],[0 1],'Color',[0.8 0.8 0.8]);
title(['Accuracy']); ylim([0 1]);  xlim([0 1]); yticks([0 0.25 0.5 0.75 1.0]); xticks([0 0.25 0.5 0.75 1.0]);

subplot(2,1,2); hold on;
scatter(reinfByDay(:,2),probeByDay(:,2),10,'filled'); plot([0 1],[0 1],'Color',[0.8 0.8 0.8]);
title(['Choice bias']); ylim([0 1.0]); yticks([0 0.25 0.5 0.75 1.0]); xlim([0 1]); xticks([0 0.25 0.5 0.75 1.0]);
end
function plotBefAft(probeCount)
%plot by day
probeByDay = fn_cell2mat(cellfun(@(x)((abs(x.probe))),probeCount{2},'UniformOutput',false),1); 
befByDay = fn_cell2mat(cellfun(@(x)((abs(x.bef))),probeCount{2},'UniformOutput',false),1); 
aftByDay = fn_cell2mat(cellfun(@(x)((abs(x.aft))),probeCount{2},'UniformOutput',false),1); 

figure;
subplot(2,1,1);
fn_boxplot({befByDay(:,1),probeByDay(:,1),aftByDay(:,1)},'paired',true,'dotLoc','none'); xticklabels({'bef','probe','aft'})
subplot(2,1,2);
fn_boxplot({befByDay(:,2),probeByDay(:,2),aftByDay(:,2)},'paired',true,'dotLoc','none'); xticklabels({'bef','probe','aft'})

[p,h,stats] = signrank(befByDay(:,1), probeByDay(:,1),'alpha',0.05/3);disp(['bef perf ' num2str(p)])
[p,h,stats] = signrank(aftByDay(:,1), probeByDay(:,1),'alpha',0.05/3);disp(['aft perf ' num2str(p)])
[p,h,stats] = signrank(befByDay(:,1), aftByDay(:,1),'alpha',0.05/3);disp(['bef aft perf ' num2str(p)])

[p,h,stats] = signrank(befByDay(:,2), probeByDay(:,2),'alpha',0.05/3);disp(['bef bias ' num2str(p)])
[p,h,stats] = signrank(aftByDay(:,2), probeByDay(:,2),'alpha',0.05/3);disp(['aft bias ' num2str(p)])
[p,h,stats] = signrank(befByDay(:,2), aftByDay(:,2),'alpha',0.05/3);disp(['bef aft bias ' num2str(p)])

fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.5);
subplot(2,1,1); hold on;
scatter(befByDay(:,1),probeByDay(:,1),10,matlabColors(1),'filled');
scatter(aftByDay(:,1),probeByDay(:,1),10,matlabColors(3),'filled'); plot([0 1],[0 1],'Color',[0.8 0.8 0.8]);
legend({'before','after'}); title('Accuracy'); ylim([0 1.0]); yticks([0 0.25 0.5 0.75 1.0]); xlim([0 1]); xticks([0 0.25 0.5 0.75 1.0]);
subplot(2,1,2); hold on;
scatter(befByDay(:,2),probeByDay(:,2),10,matlabColors(1),'filled');
scatter(aftByDay(:,2),probeByDay(:,2),10,matlabColors(3),'filled'); plot([0 1],[0 1],'Color',[0.8 0.8 0.8]);
legend({'before','after'}); title('Choice bias'); ylim([0 1.0]); yticks([0 0.25 0.5 0.75 1.0]); xlim([0 1]); xticks([0 0.25 0.5 0.75 1.0]);

fn_figureSmartDim('hSize',0.3,'widthHeightRatio',0.9); hold on;
tempBef = befByDay(:,1:2)-probeByDay(:,1:2); tempAft = aftByDay(:,1:2)-probeByDay(:,1:2);
scatter(tempBef(:,2),-tempBef(:,1),10,matlabColors(1),'filled');
scatter(tempAft(:,2),-tempAft(:,1),10,matlabColors(3),'filled'); 

lm1 = fitlm(tempBef(:,2),-tempBef(:,1)); lm2 = fitlm(tempAft(:,2),-tempAft(:,1));
tempX1 = -1:0.01:1; tempY1 = tempX1 * lm1.Coefficients{2,1} + lm1.Coefficients{1,1};
tempX2 = -1:0.01:1; tempY2 = tempX2 * lm2.Coefficients{2,1} + lm2.Coefficients{1,1};
plot(tempX1,tempY1,'Color',matlabColors(1,0.6),'LineWidth',0.5);
plot(tempX2,tempY2,'Color',matlabColors(3,0.6),'LineWidth',0.5);
legend({'before','after'}); title('Accuracy-Bias correlation'); ylim([-0.5 0.5]); yticks(-0.5:0.25:0.5); xlim([-0.8 0.8]);xticks(-0.8:0.4:0.8)

figure;
subplot(2,1,1); hold on;
p = fn_plotComparison({befByDay(:,2),probeByDay(:,2)},'compType', 'errorbarWithDot','paired',true,'dotType','random'); 
title(['Choice bias, p=' num2str(p)]); ylim([0 1]); yticks(0:0.2:1); xlim([0.6 2.4])

subplot(2,1,2); hold on;
p = fn_plotComparison({aftByDay(:,2),probeByDay(:,2)},'compType', 'errorbarWithDot','paired',true,'dotType','random'); 
title(['Choice bias, p=' num2str(p)]); ylim([0 1]); yticks(0:0.2:1); xlim([0.6 2.4])

end

