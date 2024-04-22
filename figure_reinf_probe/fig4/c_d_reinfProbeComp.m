%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);

probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;

%% 2c -- load and align probe data
% align data by minimizing performance variance
outCell = mouseMega.objFun('binProbeByTrialFromLearningOnset',{[reinfThre probeThre],probeTrialBin});
% probe data
probeData = fn_catStructField(2, outCell{1}{:});
probeIdx = cellfun(@fn_getBinMidPoint,outCell{2},'UniformOutput',false);
probeAlignPoint = cell2mat(outCell{4}); 
probeData.probeData = reshape(probeData.probeData,[size(probeData.probeData,1) 6 nMouse]);
probeData.befData = reshape(probeData.befData,[size(probeData.probeData,1) 6 nMouse]);
probeData.aftData = reshape(probeData.aftData,[size(probeData.probeData,1) 6 nMouse]);
probeData.befAftData = reshape(probeData.aftData,[size(probeData.befAftData,1) 6 nMouse]);

[probeDataRawAligned] = getProbeDataCountRaw(outCell{1},probeAlignPoint,-2:1);
[probeDataAlign,newProbeAlignPoint] = attachNan(probeData, probeAlignPoint);
% reinf data
reinfData = mouseMega.loadReinf;
reinfAlignPoint = cell2mat(outCell{3});
[reinfDataAlign,newReinfAlignPoint] = attachNan(reinfData, reinfAlignPoint);
% plot axis
reinfAxis = (1:size(reinfDataAlign.acc,1)) - newReinfAlignPoint;
tempProbeBin = probeTrialBin;
probeAxis = (tempProbeBin:tempProbeBin:tempProbeBin*size(probeDataAlign.probeAcc,1)) - tempProbeBin*newProbeAlignPoint + tempProbeBin/2;

xlimit = [-700 700]; probeSel = newProbeAlignPoint + (-2:1);

figure; subplot(2,1,1)
plotReinfProbe(reinfDataAlign.acc,probeDataAlign.probeAcc,reinfAxis,probeAxis,probeSel)
plot(xlimit,[0.5 0.5],'LineWidth',2,'Color',[0.8 0.8 0.8])
xlim(xlimit); ylim([0.5 1.0]); yticks(0.5:0.25:1.0); xticks(-600:300:600)

subplot(2,1,2)
plotReinfProbe(reinfDataAlign.bias,abs(probeDataAlign.probeBias),reinfAxis,probeAxis,probeSel)
xlim(xlimit); ylim([0 0.5]); yticks(0:0.1:0.6);xticks(-600:300:600)

%% 2c -- load and align probe data
% align data by minimizing performance variance
outCell = mouseMega.objFun('binProbeByTrialFromLearningOnset',{[reinfThre probeThre],probeTrialBin});
% probe data
probeData = fn_catStructField(2, outCell{1}{:});
probeIdx = cellfun(@fn_getBinMidPoint,outCell{2},'UniformOutput',false);
probeAlignPoint = cell2mat(outCell{4}); 
probeData.probeData = reshape(probeData.probeData,[size(probeData.probeData,1) 6 nMouse]);
probeData.befData = reshape(probeData.befData,[size(probeData.probeData,1) 6 nMouse]);
probeData.aftData = reshape(probeData.aftData,[size(probeData.probeData,1) 6 nMouse]);
probeData.befAftData = reshape(probeData.aftData,[size(probeData.befAftData,1) 6 nMouse]);

[probeDataRawAligned] = getProbeDataCountRaw(outCell{1},probeAlignPoint,-2:1);
[probeDataAlign,newProbeAlignPoint] = attachNan(probeData, probeAlignPoint);
% reinf data
reinfData = mouseMega.loadReinf;
reinfAlignPoint = cell2mat(outCell{3});
[reinfDataAlign,newReinfAlignPoint] = attachNan(reinfData, reinfAlignPoint);
% plot axis
reinfAxis = (1:size(reinfDataAlign.acc,1)) - newReinfAlignPoint;
tempProbeBin = probeTrialBin;
probeAxis = (tempProbeBin:tempProbeBin:tempProbeBin*size(probeDataAlign.probeAcc,1)) - tempProbeBin*newProbeAlignPoint + tempProbeBin/2;

xlimit = [-800 800]; probeSel = newProbeAlignPoint + (-2:1);

figure; subplot(2,1,1)
plotReinfProbe(reinfDataAlign.acc,probeDataAlign.probeAcc,reinfAxis,probeAxis,probeSel)
plot(xlimit,[0.5 0.5],'LineWidth',2,'Color',[0.8 0.8 0.8])
xlim(xlimit); ylim([0.4 1.0]); yticks(0.5:0.25:1.0); xticks(-800:400:800)

subplot(2,1,2)
plotReinfProbe(reinfDataAlign.bias,abs(probeDataAlign.probeBias),reinfAxis,probeAxis,probeSel)
xlim(xlimit); ylim([0 0.5]); yticks(0:0.1:0.6);xticks(-800:400:800)

%% 2d -- INDIVIDUAL ANIMALS
outCell = mouseMega.objFun('binProbeByDay',{}); outCell = outCell{1};
[reinfBias,reinfAcc] = cellfun(@(x)(fn_getAccBiasSmooth(x.behav.stimulus,x.behav.responseType,300)),mouseMega.mouseCell,'UniformOutput',false);
reinfData.acc = fn_cell2matFillNan(reinfAcc);
reinfData.bias = abs(fn_cell2matFillNan(reinfBias));
[reinfDataAlign,newReinfAlignPoint] = attachNan(reinfData, reinfAlignPoint);
for i = 1:nMouse
    maxidx = min([find(~isnan(reinfDataAlign.acc(:,i)),1,'last') - newReinfAlignPoint 400]);
    minidx = max([find(~isnan(reinfDataAlign.acc(:,i)),1,'first') - newReinfAlignPoint -400]);
    tempProbeIdx = outCell{i}.trial - reinfAlignPoint(i); probeSelIdx = tempProbeIdx >= minidx & tempProbeIdx<= maxidx;

    
    fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.48);
    subplot(2,1,1); hold on;  
    plot([0 0],[-5 5],'--','Color',[0.8 0.8 0.8],'LineWidth',2);
    plot(reinfAxis,(reinfDataAlign.acc(:,i)),'Color',fn_wheelColorsPT('Reinf'),'LineWidth',2);
    plot(tempProbeIdx(probeSelIdx),outCell{i}.acc(probeSelIdx),'--o','Color',fn_wheelColorsPT('Probe'),'LineWidth',2);
    xlim([minidx maxidx]); 
    if i == 10 || i == 12; ylim([0.25 1.0]); yticks([0.25 0.5 0.75 1.0]); else; ylim([0.4 1.0]); yticks([0.5 0.75 1.0]); end 
    
    subplot(2,1,2); hold on;
    plot([0 0],[-5 5],'--','Color',[0.8 0.8 0.8],'LineWidth',2);
    plot(reinfAxis,smoothdata(reinfDataAlign.bias(:,i),'gaussian',20),'Color',fn_wheelColorsPT('Reinf'),'LineWidth',2);
    plot(tempProbeIdx(probeSelIdx),abs(outCell{i}.bias(probeSelIdx)),'--o','Color',fn_wheelColorsPT('Probe'),'LineWidth',2);
    xlim([minidx maxidx]); ylim([0 1.0]); yticks([0 0.25 0.5 0.75 1.0]);
end
%% 2c -- Subsample within DAY
probeDay = mouseMega.objFun('computeProbeAlignByTrial',{'day',[reinfThre probeThre],probeTrialBin});
%% 2c -- reinf vs. probe comparison
plotSubsampleAnimal(probeDay);

%% 2d -- Combine all DAY
probeDay = mouseMega.objFun('computeProbeAlignByTrial',{'day',[reinfThre probeThre],probeTrialBin});
%% 2c -- reinf vs. probe comparison
reinfBias = fn_cell2matFillNan(cellfun(@(x)(abs(x.reinfBias)),probeDay{1},'UniformOutput',false));
probeBias = fn_cell2matFillNan(cellfun(@(x)(abs(x.probeBias)),probeDay{1},'UniformOutput',false));
reinfCri = fn_cell2matFillNan(cellfun(@(x)(abs(x.reinfCri)),probeDay{1},'UniformOutput',false));
probeCri = fn_cell2matFillNan(cellfun(@(x)(abs(x.probeCri)),probeDay{1},'UniformOutput',false));
reinfDP = fn_cell2matFillNan(cellfun(@(x)(abs(x.reinfDP)),probeDay{1},'UniformOutput',false));
probeDP = fn_cell2matFillNan(cellfun(@(x)(abs(x.probeDP)),probeDay{1},'UniformOutput',false));
reinfAcc = fn_cell2matFillNan(cellfun(@(x)(abs(x.reinfAcc)),probeDay{1},'UniformOutput',false));
probeAcc = fn_cell2matFillNan(cellfun(@(x)(abs(x.probeAcc)),probeDay{1},'UniformOutput',false));

probeAcc = fn_cell2matFillNan(cellfun(@(x)(abs(x.probeAcc)),probeDay{1},'UniformOutput',false));
probeBias = fn_cell2matFillNan(cellfun(@(x)(abs(x.probeBias)),probeDay{1},'UniformOutput',false));

figure;subplot(2,1,1);p = fn_plotComparison({reinfAcc,probeAcc},'compType', 'errorbarWithDot','paired',true);
subplot(2,1,2);p = fn_plotComparison({reinfBias,probeBias},'compType', 'errorbarWithDot','paired',true);
[p,~] = signrank(reinfBias,probeBias)
figure;subplot(2,1,1);p = fn_plotComparison({reinfDP,probeDP},'compType', 'errorbarWithDot','paired',true,'test','signrank');
subplot(2,1,2);p = fn_plotComparison({reinfCri,probeCri},'compType', 'errorbarWithDot','paired',true,'test','signrank');
%% FUNCTIONS
function plotSubsampleAnimal(probeCount)

% plot by animal
probeAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,1)'))),probeCount{2},'UniformOutput',false)); 
probeBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,2)'))),probeCount{2},'UniformOutput',false)); 
reinfAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.reinf(:,1)'))),probeCount{2},'UniformOutput',false)); 
reinfBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.reinf(:,2)'))),probeCount{2},'UniformOutput',false)); 

fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.4); 
subplot(2,1,1);hold on;
p = fn_plotComparison({nanmean(reinfAcc,2),nanmean(probeAcc,2)},'compType', 'errorbarWithDot','paired',true); 
title(['Accuracy, p=' num2str(p)]); ylim([0.5 1]);  xlim([0.6 2.4]); yticks([0.5 0.75 1.0]);

subplot(2,1,2); hold on;
p = fn_plotComparison({nanmean(reinfBias,2),nanmean(probeBias,2)},'compType', 'errorbarWithDot','paired',true); 
title(['Choice bias, p=' num2str(p)]); ylim([0 0.6]); yticks(0:0.2:0.6); xlim([0.6 2.4])

end

function [probeDataRawAligned] = getProbeDataCountRaw(probeData,alignPoint,sumIdx)

    for i = 1:length(probeData)
        tempIdx = alignPoint(i) + sumIdx; tempIdx(tempIdx<=0) = [];
        tempLen = size(probeData{i}.probeData,1); tempIdx(tempIdx>tempLen) = [];
    
        tempData.probeData(i,:) = squeeze(nansum(probeData{i}.probeData(tempIdx,:),1));
        tempData.befData(i,:) = squeeze(nansum(probeData{i}.befData(tempIdx,:),1));
        tempData.aftData(i,:) = squeeze(nansum(probeData{i}.aftData(tempIdx,:,:),1));
        tempData.befAftData(i,:) = squeeze(nansum(probeData{i}.befAftData(tempIdx,:),1));
    
    end
    
    [tempData.probeAcc,tempData.probeBias,tempData.probeDP] = fn_getAccBias(tempData.probeData); 
    [tempData.befAcc,tempData.befBias,tempData.befDP] = fn_getAccBias(tempData.befData); 
    [tempData.aftAcc,tempData.aftBias,tempData.aftDP] = fn_getAccBias(tempData.aftData); 
    [tempData.befAftAcc,tempData.befAftBias,tempData.befAftDP] = fn_getAccBias(tempData.befAftData); 
    probeDataRawAligned = tempData;
    function [acc,bias,dprime] = fn_getAccBias(mat)
        acc1 = mat(:,1)./sum(mat(:,1:2),2); acc2 = mat(:,4)./sum(mat(:,4:5),2);
        acc = (acc1+acc2)/2; bias = acc1-acc2;
        
        dprimeAccThreHigh = 0.95; dprimeAccThreLow = 0.05;
        acc1(acc1>dprimeAccThreHigh) = dprimeAccThreHigh; acc2(acc2>dprimeAccThreHigh) = dprimeAccThreHigh;
        acc1(acc1<dprimeAccThreLow) = dprimeAccThreLow; acc2(acc2<dprimeAccThreLow) = dprimeAccThreLow;
        dprime = norminv(acc1,0,1) - norminv(1-acc2,0,1);
    end
end 

function [probeData,maxAlignPoint] = attachNan(probeData, alignPoint)
    attachDim = 1; 
    maxAlignPoint = max(alignPoint);
    tempFieldNames = fieldnames(probeData);
    for i = 1:length(tempFieldNames)
        tempField = probeData.(tempFieldNames{i});
        if isnumeric(tempField)
            if length(size(tempField))<=2
                tempMatSize = size(tempField);
                tempMatSize(attachDim) = tempMatSize(attachDim) + maxAlignPoint;
                tempMat = nan(tempMatSize);
                for j = 1:length(alignPoint) % loop through the number of mouse, which is the 2nd dimension
                    startPoint = maxAlignPoint-alignPoint(j)+1;
                    tempMat(startPoint:startPoint+size(tempField,attachDim)-1,j) = tempField(:,j);
    
                end
                probeData.(tempFieldNames{i}) = tempMat;
            else
                disp(['ignore ' tempFieldNames{i} ' because >2 dimensions'])
            end
        end        
    end
end

function plotReinfProbe(reinf,probe,reinfAxis,probeAxis,probeSel)
    hold on;  
    reinfAcc_mean = nanmean(reinf,2); reinfAcc_sem = fn_nansem(reinf,2);
    probeAcc_mean = nanmean(probe(probeSel,:),2); probeAcc_sem = fn_nansem(probe(probeSel,:),2);
    plot([0 0],[-5 5],'--','Color',[0.8 0.8 0.8],'LineWidth',2);
    plot(reinfAxis,reinfAcc_mean,'Color',fn_wheelColorsPT('Reinf'),'LineWidth',2);
    f_errorbar = fn_plotFillErrorbar(reinfAxis,reinfAcc_mean,reinfAcc_sem,...
            fn_wheelColorsPT('Reinf'),'faceAlpha',0.3,'LineStyle','none');
    
    plot(probeAxis(probeSel),probeAcc_mean,'--o','Color',fn_wheelColorsPT('Probe'),'LineWidth',2);
    f_errorbar = fn_plotFillErrorbar(probeAxis(probeSel),probeAcc_mean,probeAcc_sem,...
            fn_wheelColorsPT('Probe'),'faceAlpha',0.3,'LineStyle','none'); 
end

function plotReinfProbeIndiv(reinf,probe,reinfAxis,probeAxis,probeSel)
hold on;  
plot([0 0],[-5 5],'--','Color',[0.8 0.8 0.8],'LineWidth',2);
plot(reinfAxis,reinf,'Color',fn_wheelColorsPT('Reinf'),'LineWidth',2);
plot(probeAxis(probeSel),probe(probeSel),'--o','Color',fn_wheelColorsPT('Probe'),'LineWidth',2);
end