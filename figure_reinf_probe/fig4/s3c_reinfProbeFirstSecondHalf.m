%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);

probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;

%% PLOT 1.0 -- Subsample within DAY
probeDay = mouseMega.objFun('computeProbeAlignByTrial',{'day',[reinfThre probeThre],probeTrialBin});
%% PLOT 1.1 -- reinf vs. probe comparison
plotSubsampleAnimal(probeDay);

%% FUNCTIONS
function plotSubsampleAnimal(probeCount)

% plot by animal
probeAccFirstHalf = fn_cell2matFillNan(cellfun(@(x)(abs(x.probeAccFirstHalf)),probeCount{1},'UniformOutput',false)); 
probeBiasFirstHalf = fn_cell2matFillNan(cellfun(@(x)(abs(x.probeBiasFirstHalf)),probeCount{1},'UniformOutput',false)); 
probeAccLastHalf = fn_cell2matFillNan(cellfun(@(x)(abs(x.probeAccLastHalf)),probeCount{1},'UniformOutput',false)); 
probeBiasLastHalf = fn_cell2matFillNan(cellfun(@(x)(abs(x.probeBiasLastHalf)),probeCount{1},'UniformOutput',false)); 
reinfAcc = fn_cell2matFillNan(cellfun(@(x)(abs(x.reinfAcc)),probeCount{1},'UniformOutput',false)); 
reinfBias = fn_cell2matFillNan(cellfun(@(x)(abs(x.reinfBias)),probeCount{1},'UniformOutput',false)); 

fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.4); 
subplot(3,1,1);hold on;
p = fn_plotBar({reinfAcc,probeAccLastHalf},'plotType', 'noBar','paired',true); 
title(['Accuracy, p=' num2str(p)]); ylim([0.5 1]);  xlim([0.6 2.4]); yticks([0.5 0.75 1.0]);

subplot(3,1,2); hold on;
p = fn_plotBar({reinfBias,probeBiasLastHalf},'plotType', 'noBar','paired',true); 
title(['Choice bias, p=' num2str(p)]); ylim([0 0.6]); yticks(0:0.2:0.6); xlim([0.6 2.4])

subplot(3,1,3); hold on;
p = fn_plotBar({reinfBias,probeBiasFirstHalf},'plotType', 'noBar','paired',true); 
title(['Choice bias, p=' num2str(p)]); ylim([0 0.6]); yticks(0:0.2:0.6); xlim([0.6 2.4])
subplot(3,1,3); hold on;
p = fn_plotBar({reinfBias,probeBiasFirstHalf/2+probeBiasLastHalf/2},'plotType', 'noBar','paired',true); 
title(['Choice bias, p=' num2str(p)]); ylim([0 0.6]); yticks(0:0.2:0.6); xlim([0.6 2.4])

end