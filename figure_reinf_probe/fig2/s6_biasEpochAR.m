%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT_bin40(mice);
learningCurveBin = 40;
mouseMega = wheel2AFCmega(allMouse);

%%
trialFromStartOfDay = [];
ar_bias = []; ar_unbias = [];
biasTransStart = []; unbiasTransStart = [];
for i = 1:mouseMega.nMouse
    tempAR_bias = []; tempAR_unbias = [];
    tempTrans = find(diff(allMouse{i}.biasBlock.stateFlag,1)); 
    
    for j = 1:length(tempTrans)
    if allMouse{i}.biasBlock.stateFlag(tempTrans(j)+1)~=0
        biasTransStart = [biasTransStart allMouse{i}.behav.masterTrialnumSession(tempTrans(j)) ];
    else
        unbiasTransStart = [unbiasTransStart allMouse{i}.behav.trialnumSession(tempTrans(j)) ];

    end

    end 


    for j = 1:length(allMouse{i}.biasBlock.blockL.start)
        tempTrial = allMouse{i}.biasBlock.blockL.start(j);
        tempDay = allMouse{i}.behav.day(tempTrial); 
        tempTrialDayStart = tempTrial - find(allMouse{i}.behav.day==tempDay,1); 
        trialFromStartOfDay = [trialFromStartOfDay tempTrialDayStart];
        nTrials = 30;
        try
            tempBiasBlockFlag1 = [allMouse{i}.biasBlock.blockL.start(j)+nTrials allMouse{i}.biasBlock.blockL.end(j)-nTrials];
            tempBiasBlockFlag2 = [allMouse{i}.biasBlock.blockL.start(j)-nTrials allMouse{i}.biasBlock.blockL.end(j)+nTrials];
            tempFlag2Sel = allMouse{i}.biasBlock.stateFlag(tempBiasBlockFlag2)==0;
            tempAR1 = allMouse{i}.behav.actionRate(tempBiasBlockFlag1)';
            tempAR2 = allMouse{i}.behav.actionRate(tempBiasBlockFlag2)';
            tempAR2 = tempAR2(tempFlag2Sel);
            
        catch
            disp('a1');
            tempAR1 = []; tempAR2 = [];
        end
        tempAR_bias = [tempAR_bias tempAR1];
        tempAR_unbias = [tempAR_unbias tempAR2];
    end

    for j = 1:length(allMouse{i}.biasBlock.blockR.start)
        tempTrial = allMouse{i}.biasBlock.blockR.start(j);
        tempDay = allMouse{i}.behav.day(tempTrial); 
        tempTrialDayStart = tempTrial - find(allMouse{i}.behav.day==tempDay,1); 
        trialFromStartOfDay = [trialFromStartOfDay tempTrialDayStart];
    

        try
            tempBiasBlockFlag1 = [allMouse{i}.biasBlock.blockR.start(j)+nTrials allMouse{i}.biasBlock.blockR.end(j)-nTrials];
            tempBiasBlockFlag2 = [allMouse{i}.biasBlock.blockR.start(j)-nTrials allMouse{i}.biasBlock.blockR.end(j)+nTrials];
            tempFlag2Sel = allMouse{i}.biasBlock.stateFlag(tempBiasBlockFlag2)==0;
            tempAR1 = allMouse{i}.behav.actionRate(tempBiasBlockFlag1)';
            tempAR2 = allMouse{i}.behav.actionRate(tempBiasBlockFlag2)';
            tempAR2 = tempAR2(tempFlag2Sel);
            
        catch
            disp('a2');
            tempAR1 = []; tempAR2 = [];
        end
        tempAR_bias = [tempAR_bias tempAR1];
        tempAR_unbias = [tempAR_unbias tempAR2];
    end
    ar_bias = [ar_bias nanmean(tempAR_bias)];
    ar_unbias = [ar_unbias nanmean(tempAR_unbias)];

end 
figure; 
subplot(2,1,1);fn_plotHistLine(biasTransStart,'histCountArgIn',{0:25:300,'Normalization','count'}); 
xlim([-100 400]); ylim([0 20]); title('transition into biased epoch')
subplot(2,1,2);fn_plotHistLine(unbiasTransStart,'histCountArgIn',{0:25:300,'Normalization','count'}); 
xlim([-100 400]); ylim([0 20]); title('transition into unbiased epoch')
figure; fn_plotComparison({ar_bias,ar_unbias},'paired', true); ylim([0.5 1]); 
%%
for i = 1:mouseMega.nMouse
    figure; windowSize = 400;
    tempProbeIdx = mouseMega.mouseCell{i}.probe.onIdx;
    for j = 1:20
        subplot(4,5,j); hold on; 
        tempStart = (j-1)*windowSize+1:1:j*windowSize; 
        if tempStart(end) < size(mouseMega.mouseCell{i}.behav,1)
            probeSelIdx = tempProbeIdx > tempStart(1) & tempProbeIdx < tempStart(end);
    
            scatter(mouseMega.mouseCell{i}.behav.acc(tempStart),...
                mouseMega.mouseCell{i}.behav.bias(tempStart),10,matlabColors(1),'filled');
            plot(mouseMega.mouseCell{i}.behav.acc(tempStart),...
                mouseMega.mouseCell{i}.behav.bias(tempStart),'Color',matlabColors(1));
            if any(probeSelIdx)
                tempProbe = sum(mouseMega.mouseCell{i}.probe.probeData(probeSelIdx,:),1);
                tempReinf = sum(mouseMega.mouseCell{i}.probe.befData(probeSelIdx,:) + mouseMega.mouseCell{i}.probe.aftData(probeSelIdx,:),1);
                [pBias,pAcc] = fn_getAccBiasByCount(tempProbe);
                [rBias,rAcc] = fn_getAccBiasByCount(tempReinf);
                scatter([pAcc],[pBias],20,matlabColors(2),'filled');
                scatter([rAcc],[rBias],20,matlabColors(6),'filled');
            end
            xlim([0.2 1]); ylim([-1 1])
        end
    end 

end

%%

for i = 1:mouseMega.nMouse
    figure; windowSize = 200;
    tempProbeIdx = mouseMega.mouseCell{i}.probe.onIdx;
    for j = 1:12
        subplot(3,4,j); hold on; 
        tempStart = (j-1)*windowSize+1:2:j*windowSize; 
        
        probeSelIdx = tempProbeIdx > tempStart(1) & tempProbeIdx < tempStart(end);
        try
            scatter(mouseMega.mouseCell{i}.behav.acc1(tempStart),...
                mouseMega.mouseCell{i}.behav.acc2(tempStart),10,matlabColors(1),'filled');
            plot(mouseMega.mouseCell{i}.behav.acc1(tempStart),...
                mouseMega.mouseCell{i}.behav.acc2(tempStart),'Color',matlabColors(1));
            if ~isempty(probeSelIdx)
                tempProbe = sum(mouseMega.mouseCell{i}.probe.probeData(probeSelIdx,:));
                tempReinf = sum(mouseMega.mouseCell{i}.probe.befData(probeSelIdx,:) + mouseMega.mouseCell{i}.probe.aftData(probeSelIdx,:));
                [pBias,pAcc,~,pAcc1,pAcc2] = fn_getAccBiasByCount(tempProbe);
                [rBias,rAcc,~,rAcc1,rAcc2] = fn_getAccBiasByCount(tempReinf);
                scatter([pAcc1],[pAcc2],20,matlabColors(2),'filled');
                scatter([rAcc1],[rAcc2],20,matlabColors(6),'filled');
            end
        catch
        end
        xlim([0 1]); ylim([0 1])
    end 

end