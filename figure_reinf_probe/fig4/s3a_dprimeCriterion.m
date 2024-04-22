clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);

probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;

%% Subsample within DAY
probeDay = mouseMega.objFun('computeProbeAlignByTrial',{'day',[reinfThre probeThre],probeTrialBin});
%% supplement -- reinf vs. probe comparison
plotSubsampleAnimal(probeDay);
%% supplement -- reinf vs. probe comparison
plotScatter(probeDay);
%% Subsample within BEF and AFT
probeDay = mouseMega.objFun('computeProbeAlignByTrial',{'befAft',[reinfThre probeThre],probeTrialBin});
%% supplement -- reinf vs. probe comparison
plotSubsampleAnimal(probeDay);
%% FUNCTIONS
function plotSubsampleAnimal(probeCount)
% plot by animal
probeAcc = fn_cell2matFillNan(cellfun(@(x)(((x.probe(:,3)'))),probeCount{2},'UniformOutput',false)); 
probeBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,4)'))),probeCount{2},'UniformOutput',false)); 
reinfAcc = fn_cell2matFillNan(cellfun(@(x)(((x.reinf(:,3)'))),probeCount{2},'UniformOutput',false)); 
reinfBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.reinf(:,4)'))),probeCount{2},'UniformOutput',false)); 

fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.4); 
p = fn_plotComparison({reinfBias(:),probeBias(:)},'compType', 'errorbarWithDot','paired',true); 

fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.4); 
subplot(2,1,1);hold on;
p = fn_plotComparison({nanmean(reinfAcc,2),nanmean(probeAcc,2)},'compType', 'errorbarWithDot','paired',true); 
title(['Accuracy, p=' num2str(p)]); ylim([0 3]);  xlim([0.6 2.4]); yticks(0:1:3);

subplot(2,1,2); hold on;
p = fn_plotComparison({nanmean(reinfBias,2),nanmean(probeBias,2)},'compType', 'errorbarWithDot','paired',true); 
title(['Choice bias, p=' num2str(p)]); ylim([0 1]); xlim([0.6 2.4]); yticks(0:0.5:1.0); 
end
function plotScatter(probeCount)

% plot by animal
probeAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,3)'))),probeCount{2},'UniformOutput',false)); 
probeBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,4)'))),probeCount{2},'UniformOutput',false)); 
reinfAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.reinf(:,3)'))),probeCount{2},'UniformOutput',false)); 
reinfBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.reinf(:,4)'))),probeCount{2},'UniformOutput',false)); 

colors = parula; tempIdx = round(linspace(1, size(colors,1),size(reinfBias,1)+1));
colors = colors(tempIdx(1:end-1),:);
tempColor = []; tempColor(:,1,:) = colors; tempColor = repmat(tempColor,1,size(probeAcc,2),1);
tempColor = reshape(tempColor,size(tempColor,1)* size(tempColor,2),[]);

fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.5); 
subplot(2,1,1); tempB = reinfBias(:)-probeBias(:); tempA = probeAcc(:)-reinfAcc(:); 
scatter(reinfBias(:)-probeBias(:),probeAcc(:)-reinfAcc(:),30,tempColor,'filled'); hold on;
lm = fitlm(tempB,tempA);
tempX = -1:0.01:1; tempY = tempX * lm.Coefficients{2,1} + lm.Coefficients{1,1};
plot(tempX,tempY,'Color',[0.6 0.6 0.6],'LineWidth',0.5);
xlim([-1 1]); xticks(-1:0.5:1); ylim([-1.5 2]);yticks(-1.5:0.5:2);
title(['corr = ' num2str(corr(tempB,tempA,'rows','complete')) ', p = ' num2str(lm.Coefficients.pValue(2))])

subplot(2,1,2); tempB = reinfBias(:); tempA = reinfBias(:)-probeBias(:);
scatter(tempB, tempA,30,tempColor,'filled'); hold on;
lm = fitlm(tempB,tempA);
tempX = 0:0.01:1.5; tempY = tempX * lm.Coefficients{2,1} + lm.Coefficients{1,1};
plot(tempX,tempY,'Color',[0.6 0.6 0.6],'LineWidth',0.5);
xlim([0 1.5]); xticks(0:0.5:1.5); ylim([-1 1]);yticks(-1:0.5:1)
title(['corr = ' num2str(corr(tempB,tempA,'rows','complete')) ', p = ' num2str(lm.Coefficients.pValue(2))])
end