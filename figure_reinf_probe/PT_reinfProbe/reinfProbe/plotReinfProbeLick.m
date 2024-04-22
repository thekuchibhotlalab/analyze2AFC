%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};

%mice ={'zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);


probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
%% check lick dection method
for i = 1:length(allMouse)
    nDay = unique(allMouse{i}.behav.day);
    for j = 1:length(nDay)
        nLick{i}(j) = sum(allMouse{i}.behav.nlick(allMouse{i}.behav.day==nDay(j)));
    end

end



%% 
lickData = mouseMega.objFun('computeLickAlignByTrial',{'day',[reinfThre probeThre],probeTrialBin});
lickReinf = []; lickProbe = []; lickTimeReinf = []; lickTimeProbe = [];
lickReinf(:,1) = fn_cell2mat(cellfun(@(x)(nanmean(x{1},2)),lickData{1},'UniformOutput',false),1); 
lickReinf(:,2) = fn_cell2mat(cellfun(@(x)(nanmean(x{2},2)),lickData{1},'UniformOutput',false),1); 

lickProbe(:,1) = fn_cell2matFillNan(cellfun(@(x)(nanmean(x{1},2)),lickData{2},'UniformOutput',false)); 
lickProbe(:,2) = fn_cell2matFillNan(cellfun(@(x)(nanmean(x{2},2)),lickData{2},'UniformOutput',false)); 


lickTimeReinf(:,1) = fn_cell2mat(cellfun(@(x)(nanmean(x{1},2)),lickData{3},'UniformOutput',false),1); 
lickTimeReinf(:,2) = fn_cell2mat(cellfun(@(x)(nanmean(x{2},2)),lickData{3},'UniformOutput',false),1); 

lickTimeProbe(:,1) = fn_cell2matFillNan(cellfun(@(x)(nanmean(x{1},2)),lickData{4},'UniformOutput',false)); 
lickTimeProbe(:,2) = fn_cell2matFillNan(cellfun(@(x)(nanmean(x{2},2)),lickData{4},'UniformOutput',false)); 

%%

selMouse = [2 3 5 12 13]; figure;
fn_plotBar({lickReinf(selMouse,1),lickReinf(selMouse,2),lickProbe(selMouse,1),lickProbe(selMouse,2)});

figure; subplot(4,1,1)
lickRasterReinf = fn_cell2matFillNan(lickData{5}{2}{1}); 
xaxis = repmat(1:size(lickRasterReinf,1),size(lickRasterReinf,2),1)';
scatter(lickRasterReinf(:),xaxis(:),5,[0.8 0.8 0.8],'filled'); xlim([-2 4])

subplot(4,1,2)
lickRasterReinf = fn_cell2matFillNan(lickData{5}{2}{2}); 
xaxis = repmat(1:size(lickRasterReinf,1),size(lickRasterReinf,2),1)';
scatter(lickRasterReinf(:),xaxis(:),5,[0.8 0.8 0.8],'filled'); xlim([-2 4])

subplot(4,1,3)
lickRasterReinf = fn_cell2matFillNan(lickData{6}{2}{1}); 
xaxis = repmat(1:size(lickRasterReinf,1),size(lickRasterReinf,2),1)';
scatter(lickRasterReinf(:),xaxis(:),5,[0.8 0.8 0.8],'filled'); xlim([-2 4])

subplot(4,1,4)
lickRasterReinf = fn_cell2matFillNan(lickData{6}{2}{2}); 
xaxis = repmat(1:size(lickRasterReinf,1),size(lickRasterReinf,2),1)';
scatter(lickRasterReinf(:),xaxis(:),5,[0.8 0.8 0.8],'filled'); xlim([-2 4])

