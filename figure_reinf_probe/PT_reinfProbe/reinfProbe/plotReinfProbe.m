function plotReinfProbe(mouseMega)

%[acc,bias,probeIdx] = mouseMega.loadAccBias;
reinf = mouseMega.loadReinf;
probeBin = mouseMega.loadProbeBin;

nMouse = size(reinf.acc,2);

trialLim = 2100; trialBin = 0:300:trialLim;
plotLim = 2000;
trialBinPlot = (trialBin(1:end-1) + trialBin(2:end))/2;
binProbe = nan(length(trialBin)-1,nMouse);
binProbeBias = nan(length(trialBin)-1,nMouse);
for i = 1:nMouse
    trialNum = probeBin.probeIdx(:,i);
    trialLimFlag = trialNum>trialLim | isnan(trialNum); trialNum(trialLimFlag) = [];
    
    binIdx = sum(trialNum>trialBin,2);
    binProbe(binIdx,i) = probeBin.probeAcc(~trialLimFlag,i);
    binProbeBias(binIdx,i) = probeBin.probeBias(~trialLimFlag,i); 
end
%% ----------------------------PLOT 1 PLOT REINF-PROBE -------------------------------
plotSEM(reinf.acc, binProbe, plotLim, trialBinPlot);ylabel('Accuracy');
plot([1 trialLim],[0.5 0.5],'Color',[0.6 0.6 0.6],'LineWidth',2);
ylim([0.4 1.0]); yticks([0.5 0.75 1])
try
    modelAcc = mouseMega.getProp('behav','field','modelAcc','matFlag',true);
    modelBias = mouseMega.getProp('behav','field','modelBias','matFlag',true);
    plotSEM(reinf.acc, binProbe, plotLim, trialBinPlot);ylabel('Accuracy');
    plot([1 trialLim],[0.5 0.5],'Color',[0.6 0.6 0.6],'LineWidth',2);
    ylim([0.4 1.0]); yticks([0.5 0.75 1])
    f_errorbar = fn_plotFillErrorbar(1:size(modelAcc,1),(nanmean(modelAcc,2))',...
            (nanstd(modelAcc,0,2)./sqrt(size(modelAcc,2)))',...
            matlabColors(4),'faceAlpha',0.3,'LineStyle','none');
    plot(nanmean(modelAcc,2),'Color', matlabColors(4), 'LineWidth',2);
catch
end

plotIndividual(reinf.acc,binProbe,plotLim,trialBinPlot);ylabel('Accuracy');
plot([1 trialLim],[0.5 0.5],'Color',[0.6 0.6 0.6],'LineWidth',2);
ylim([0.4 1.0]); yticks([0.5 0.75 1])

plotIndividualReinf(reinf.acc,plotLim);ylabel('Accuracy');
plot([1 trialLim],[0.5 0.5],'Color',[0.6 0.6 0.6],'LineWidth',2);

plotSEM(reinf.bias, abs(binProbeBias), plotLim, trialBinPlot); ylabel('Bias');
ylim([0 0.7]); yticks(0:0.2:0.6)
plotIndividual(reinf.bias, abs(binProbeBias), plotLim, trialBinPlot); ylabel('Bias');
plotIndividualReinf(reinf.bias,plotLim); ylabel('Bias');

%% ----------------------------PLOT 2 BAR REINF BEF AFT PROBE -------------------------------

blockSelFlag = selectProbeBlock(reinf,probeBin,probeBin.probeIdx); probeSelFlag = blockSelFlag;
for i = 1:size(blockSelFlag,2); tempidx = find(blockSelFlag(:,i) == 2,1); probeSelFlag(tempidx:end,i) = 0; end
disp('Number of probe blocks per animal'); disp(sum(probeSelFlag,1))

probeSelFlag = logical(probeSelFlag);
reinfAcc = (reshape(probeBin.befAcc(probeSelFlag),1,[]) + reshape(probeBin.aftAcc(probeSelFlag),1,[]))/2;
probeAcc = reshape(probeBin.probeAcc(probeSelFlag),1,[]);

modelBiasBef = probeBin.befModelBias; %mouseMega.getProp('probe','field','befModelBias','matFlag',true);
modelBiasAft = probeBin.aftModelBias; %mouseMega.getProp('probe','field','aftModelBias','matFlag',true);
reinfModelBias = (abs(reshape(modelBiasBef(probeSelFlag),1,[])) + abs(reshape(modelBiasAft(probeSelFlag),1,[])))/2;

reinfBias = abs(reshape(probeBin.befBias(probeSelFlag),1,[]) + abs(reshape(probeBin.aftBias(probeSelFlag),1,[])))/2;
probeBias = abs(reshape(probeBin.probeBias(probeSelFlag),1,[])) ;

figure; scatter(reinfModelBias,probeAcc-reinfAcc,40,'filled');xlim([0 1]);ylim([-0.3 0.5])
mdl = fitlm(reinfModelBias,probeAcc-reinfAcc);
xlabel('Model Bias'); ylabel('P-R Acc Diff')

figure; scatter(reinfModelBias-probeBias,probeAcc-reinfAcc,40,'filled');xlim([-1 1]);ylim([-0.3 0.5])
mdl = fitlm(reinfModelBias-probeBias,probeAcc-reinfAcc);
xlabel('Model Bias - Probe Bias'); ylabel('P-R Acc Diff')

figure; scatter(reinfBias-probeBias,probeAcc-reinfAcc,40,'filled');xlim([-1 1]);ylim([-0.3 0.5])
mdl = fitlm(reinfBias-probeBias,probeAcc-reinfAcc);
xlabel('reinf Bias - Probe Bias'); ylabel('P-R Acc Diff')


figure; subplot(2,1,1); hold on;
plotBar({probeBin.befAcc,probeBin.probeAcc,probeBin.aftAcc},probeSelFlag==1,'left',false,false); ylabel('Accuracy');
subplot(2,1,2); hold on;
plotBar({probeBin.befBias,probeBin.probeBias,probeBin.aftBias},probeSelFlag==1,'right',true,false); ylabel('Bias');

figure; subplot(2,1,1); hold on;
plotBar({probeBin.befAcc,probeBin.probeAcc,probeBin.aftAcc},probeSelFlag==1,'left',false,true); ylabel('Accuracy');
ylim([0.4 1])
subplot(2,1,2); hold on;
plotBar({probeBin.befBias,probeBin.probeBias,probeBin.aftBias},probeSelFlag==1,'right',true,true); ylabel('Bias');
ylim([0 0.7])

figure; plotScatterBias(probeBin,probeSelFlag==1,true);
xlabel('Reinf');ylabel('Probe'); title('Bias')

highBiasFlag = (abs(probeBin.befBias)+abs(probeBin.aftBias))/2 > 0.3;
figure; subplot(2,1,1); hold on;
plotBar({probeBin.befAcc,probeBin.probeAcc,probeBin.aftAcc},probeSelFlag==1 & highBiasFlag,'left',false,true); ylabel('Accuracy');
subplot(2,1,2); hold on;
plotBar({probeBin.befBias,probeBin.probeBias,probeBin.aftBias},probeSelFlag==1 & highBiasFlag,'right',true,true); ylabel('Bias');

figure; subplot(2,1,1); hold on;
plotBar({probeBin.befAcc,probeBin.probeAcc,probeBin.aftAcc},probeSelFlag==1 & (~highBiasFlag),'left',false,true); ylabel('Accuracy');
subplot(2,1,2); hold on;
plotBar({probeBin.befBias,probeBin.probeBias,probeBin.aftBias},probeSelFlag==1 & (~highBiasFlag),'right',true,true); ylabel('Bias');
%% plot 3 -- plot reinf vs. model diff

sampInterval = 20;
modelSelIdx = probeBin.probeIdx; modelSelIdx(~probeSelFlag) = nan;
modelAccSel_A_all = [];modelAccSel_H_all = []; accSel_all = []; modelBiasSel_all = []; biasSel_all = [];
tempStartAll = []; tempEndAll = [];
modelAccSel = nan(1,size(modelSelIdx,2));accSel = nan(1,size(modelSelIdx,2));

modelAcc_removeA = mouseMega.getProp('behav','field','modelPred_removeAAcc','matFlag',true);
modelAcc_removeH = mouseMega.getProp('behav','field','modelPred_removeHAcc','matFlag',true);
modelAcc_removeAll = mouseMega.getProp('behav','field','modelPred_removeAllAcc','matFlag',true);

% for removal of action, bias or remove all
modelAccSel_A = nan(1,size(modelSelIdx,2));
modelAccSel_H = nan(1,size(modelSelIdx,2));
modelAccSel_All = nan(1,size(modelSelIdx,2));
for i = 1:size(modelSelIdx,2)
    tempStart = round(min(modelSelIdx(:,i))-150); tempStartAll(i) = tempStart;
    tempEnd = round(max(modelSelIdx(:,i))+150);  tempEndAll(i) = tempEnd;
    %tempStart = find(smoothdata(reinf.acc(:,i),'movmean',100)>0.7,1);
    %tempEnd = find(smoothdata(reinf.acc(:,i),'movmean',100)>0.85,1);
    %if isempty(tempEnd); tempEnd = find(isnan(reinf.acc(:,i)),1)-1; end 
    if ~isnan(tempStart) && ~isnan(tempEnd)
        modelAccSel_A(i) = nanmean(modelAcc_removeA(tempStart:sampInterval:tempEnd,i),1);
        modelAccSel_H(i) = nanmean(modelAcc_removeH(tempStart:sampInterval:tempEnd,i),1);
        modelAccSel_All(i) = nanmean(modelAcc_removeAll(tempStart:sampInterval:tempEnd,i),1);

        modelAccSel(i) = nanmean(modelAcc(tempStart:sampInterval:tempEnd,i),1);
        accSel(i) = nanmean(reinf.acc(tempStart:sampInterval:tempEnd,i),1);
        
        modelAccSel_A_all = cat(1,modelAccSel_A_all,modelAcc_removeA(tempStart:sampInterval:tempEnd,i));
        modelAccSel_H_all = cat(1,modelAccSel_H_all,modelAcc_removeH(tempStart:sampInterval:tempEnd,i));
        
        accSel_all = cat(1,accSel_all,reinf.acc(tempStart:sampInterval:tempEnd,i));
        biasSel_all = cat(1,biasSel_all,bias.reinf(tempStart:sampInterval:tempEnd,i));
        modelBiasSel_all = cat(1,modelBiasSel_all,modelBias(tempStart:sampInterval:tempEnd,i));
    end
end
figure;
bar([nanmean(accSel) nanmean(modelAccSel) nanmean(modelAccSel_A) nanmean(modelAccSel_H) nanmean(modelAccSel_All)] ...
    ,'EdgeColor',[0 0 0],'FaceColor','None'); hold on;
ylim([0.4 1])


disp(['start trial:' num2str(nanmean(tempStartAll))]); disp(['end trial:' num2str(nanmean(tempEndAll))])
figure;
[h,p1] = ttest(accSel,modelAccSel_A);
[h,p2] = ttest(accSel,modelAccSel_H);
[h,p3] = ttest(modelAccSel_A,modelAccSel_H);
bar([nanmean(accSel) nanmean(modelAccSel_A) nanmean(modelAccSel_H)] ,'EdgeColor',[0 0 0],'FaceColor','None'); hold on;
for i = 1:length(modelAccSel_A)
    plot([1 2 3],[accSel(i) modelAccSel_A(i) modelAccSel_H(i)],'Color',[0.6 0.6 0.6],'Marker','.','MarkerSize',15,...
        'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',0.5);
end
legend({['pBef = ' num2str(p1,'%.2e')], ['pBef = ' num2str(p2,'%.2e') ], ['pBef = ' num2str(p3,'%.2e') ]},'Location','Best')
xticks([1 2 3]); xticklabels({'behavior','remove A', 'removeo H'}); xlim([0 4]);
ylim([0.4 1])

figure; subplot(1,2,1)
scatter(abs(modelBiasSel_all(1:end)),(modelAccSel_A_all(1:end)-accSel_all(1:end)),15,'filled')
xlabel('bias'); ylabel('Model Increase in performance')
subplot(1,2,2)
scatter(abs(biasSel_all(1:end)),(modelAccSel_A_all(1:end)-accSel_all(1:end)),15,'filled')
xlabel('bias'); ylabel('Model Increase in performance')

figure;hold on; tempx = abs(biasSel_all(1:end)); 
tempy1 = (modelAccSel_A_all(1:end)-accSel_all(1:end));  tempy2 = (modelAccSel_H_all(1:end)-accSel_all(1:end));

scatter(tempx,tempy2,10,[0.4 0.4 0.4],'filled','MarkerFaceAlpha',0.5); hold on;
scatter(tempx,tempy1,10,matlabColors(2,0.9),'filled','MarkerFaceAlpha',0.5);

lm = fitlm(tempx , tempy1); 
plot([0 1],[0 1]*lm.Coefficients{2,1}+lm.Coefficients{1,1},'LineWidth',2,'Color',matlabColors(2,0.9))
lm = fitlm(tempx , tempy2); 
plot([0 1],[0 1]*lm.Coefficients{2,1}+lm.Coefficients{1,1},'LineWidth',2,'Color',[0.4 0.4 0.4])

xlabel('bias'); ylabel('Model Increase in performance');
xlim([0 1.0]); xticks([0 0.5 1.0])
ylim([-0.1 0.25]); yticks([-0.1 0 0.1 0.2])
disp(['p = ' num2str(lm.Coefficients{2,4})])
end

function plotSEM(reinf, binProbe, trialLim, trialBinPlot)
    [meanProbe,semProbe,nanflag] = getBinProbe(binProbe);

    fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
    f_errorbar = fn_plotFillErrorbar(trialBinPlot(~nanflag),meanProbe(~nanflag),semProbe(~nanflag),...
        fn_wheelColorsPT('Probe'),'faceAlpha',0.3,'LineStyle','none');
    f_errorbar = fn_plotFillErrorbar(1:size(reinf,1),(nanmean(reinf,2))',...
        (nanstd(reinf,0,2)./sqrt(size(reinf,2)))',...
        fn_wheelColorsPT('Reinf'),'faceAlpha',0.3,'LineStyle','none');
    plot(nanmean(reinf,2),'Color', fn_wheelColorsPT('Reinf'), 'LineWidth',2);
    plot(trialBinPlot(~nanflag),meanProbe(~nanflag),'--o','Color',fn_wheelColorsPT('Probe'),'LineWidth',2)
    xlim([1 trialLim]); xlabel('Trials');
    
end

function plotIndividual(reinf,binProbe,trialLim,trialBinPlot)
    [meanProbe,semProbe,nanflag] = getBinProbe(binProbe);
    fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
    % PLOT REINF
    for j = 1:size(reinf,2)
        plot(reinf(:,j),'Color',  fn_wheelColorsPT('Reinf',0.3), 'LineWidth',1);
    end
    for j = 1:size(binProbe,2)
        plot(trialBinPlot,binProbe(:,j),'--o','Color',  fn_wheelColorsPT('Probe',0.3), 'LineWidth',1)
    end
    plot(nanmean(reinf,2),'Color',  fn_wheelColorsPT('Reinf'), 'LineWidth',2);
    % PLOT PROBE AGAIN
    plot(trialBinPlot(~nanflag),meanProbe(~nanflag),'--o','Color',fn_wheelColorsPT('Probe'),'LineWidth',2)
    %plot([1 trialLim],[0.5 0.5],'Color',[0.6 0.6 0.6],'LineWidth',2);
    xlim([1 trialLim]);  xlabel('Trials');
end

function plotIndividualReinf(reinf,trialLim)
    fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
    % PLOT REINF
    for j = 1:size(reinf,2)
        plot(reinf(:,j),'Color', fn_wheelColorsPT('Reinf',0.3), 'LineWidth',1);
    end
    plot(nanmean(reinf,2),'Color', fn_wheelColorsPT('Reinf'), 'LineWidth',2);
    xlim([1 trialLim]); xlabel('Trials');
end

function plotBar(S,blockSelFlag,tail,absFlag,combineAni)
    
    %------------------------PLOT BAR FOR EACH ANIMAL---------------------------
    if combineAni
        tempProbeFlat = S{2}; tempProbeFlat(~blockSelFlag) = nan; 
        tempBefFlat = S{1}; tempBefFlat(~blockSelFlag) = nan; 
        tempAftFlat = S{3}; tempAftFlat(~ blockSelFlag) = nan; 
        if absFlag; tempBefFlat = abs(tempBefFlat); tempAftFlat = abs(tempAftFlat); tempProbeFlat = abs(tempProbeFlat); end
        tempProbeFlat = nanmean(tempProbeFlat,1); tempBefFlat = nanmean(tempBefFlat,1); tempAftFlat = nanmean(tempAftFlat,1);
    %------------------------PLOT BAR FOR EACH SESSION---------------------------
    else 
        tempBefFlat = reshape(S{1}(blockSelFlag),1,[]);
        tempAftFlat = reshape(S{3}(blockSelFlag),1,[]);
        tempProbeFlat = reshape(S{2}(blockSelFlag),1,[]);
        if absFlag; tempBefFlat = abs(tempBefFlat); tempAftFlat = abs(tempAftFlat); tempProbeFlat = abs(tempProbeFlat); end
    end
    
    nanFlag = isnan(tempProbeFlat); 
    tempBefFlat(nanFlag) = []; tempAftFlat(nanFlag) = []; tempProbeFlat(nanFlag) = [];
    [pBef,hBef] = signrank(tempBefFlat,tempProbeFlat,'tail',tail);
    [pAft,hAft] = signrank(tempAftFlat,tempProbeFlat,'tail',tail);
    bar([nanmean(tempBefFlat) nanmean(tempProbeFlat) nanmean(tempAftFlat)] ,'EdgeColor',[0 0 0],'FaceColor','None'); hold on;
    for i = 1:length(tempBefFlat)
        f = plot([1 2 3],[tempBefFlat(i) tempProbeFlat(i) tempAftFlat(i)],'Color',[0.6 0.6 0.6],'Marker','.','MarkerSize',15,...
            'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',0.5);
    end
    legend(f,['pBef = ' num2str(pBef,'%.2e') newline 'pAft = ' num2str(pAft,'%.2e')],'Location','Best')
    xticks([1 2 3]); xticklabels({'Bef','Probe','Aft'}); xlim([0 4]);
    
end

function plotScatterBias(S,blockSelFlag,absFlag)
    tempBefFlat = reshape(S.befBias(blockSelFlag),1,[]);
    tempAftFlat = reshape(S.aftBias(blockSelFlag),1,[]);
    tempProbeFlat = reshape(S.probeBias(blockSelFlag),1,[]);
    if absFlag; tempBefFlat = abs(tempBefFlat); tempAftFlat = abs(tempAftFlat); tempProbeFlat = abs(tempProbeFlat); end
    nanFlag = isnan(tempProbeFlat); 
    tempBefFlat(nanFlag) = []; tempAftFlat(nanFlag) = []; tempProbeFlat(nanFlag) = [];
    plot([0 1],[0 1],'LineWidth',2,'Color',[0.8 0.8 0.8]);hold on;
    scatter(tempBefFlat,tempProbeFlat,50,matlabColors(3,0.7),'filled','^'); 
    scatter(tempAftFlat,tempProbeFlat,50,matlabColors(4,0.7),'filled','v');
    xlim([0 1]); ylim([0 1]); xticks([0 0.5 1]); yticks([0 0.5 1])
end

function [meanProbe,semProbe,nanflag] = getBinProbe(binProbe)
    nData = sum(~isnan(binProbe),2); nData(nData==0) = nan;
    meanProbe = nanmean(binProbe,2); nanflag = isnan(meanProbe);
    semProbe = nanstd(binProbe,0,2)./sqrt(nData);
end

function blockSelFlag = selectProbeBlock(reinf,probeBin,probeIdx)
    nTrialBefAft = 100;
    lowAccThre = 0.65; maxAccThre = 0.85;
    lowAccThreProbe = 0.65;
    probeIdx = round(probeIdx);
    blockSelFlag = zeros(size(probeIdx));
    for i = 1:size(probeIdx,2) % loop through all animals
        maxAcc = max(smoothdata(reinf.acc(:,i),'movmean',200));
        accThre =  maxAccThre * maxAcc; if accThre<0.75; accThre = 0.75; end
        for j = 1:size(probeIdx,1) % loop through all days 
            if isnan(probeIdx(j,i)); tempAcc = 0;
            else
                tempStart = probeIdx(j,i)-nTrialBefAft; tempEnd = probeIdx(j,i)+nTrialBefAft;
                if tempStart > size(reinf.acc,1); tempStart = size(reinf.acc,1); end
                if tempEnd > size(reinf.acc,1); tempEnd = size(reinf.acc,1); end
                tempAcc = nanmean(reinf.acc(tempStart:tempEnd,i));
            end
            
            if tempAcc > lowAccThre || probeBin.probeAcc(j,i) > lowAccThreProbe
                if  tempAcc > accThre; blockSelFlag(j,i) = 2;
                %if  tempAcc > maxAccThre; blockSelFlag(j,i) = 2;
                else; blockSelFlag(j,i) = 1;
                end
            end

        end
    end



end