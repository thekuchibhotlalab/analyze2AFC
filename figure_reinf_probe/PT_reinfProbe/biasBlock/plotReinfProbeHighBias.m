function plotReinfProbeHighBias(mouseMega)

stateFlag = mouseMega.getProp('biasBlock','field','stateFlag');
probeFlag = mouseMega.getProp('behav','field','probe');
befFlag = mouseMega.getProp('behav','field','reinfBef');
aftFlag = mouseMega.getProp('behav','field','reinfAft');
stim = mouseMega.getProp('behav','field','stimulus');
responseType = mouseMega.getProp('behav','field','responseType');

[ probeAcc, probeBias ]= cellfun(@(x,y,z,a)getProbe(x,y,z,a),stim,responseType,stateFlag,probeFlag,'UniformOutput', true);
[ befAcc,befBias ]= cellfun(@(x,y,z,a)getProbe(x,y,z,a),stim,responseType,stateFlag,befFlag,'UniformOutput', true);
[aftAcc,aftBias ]= cellfun(@(x,y,z,a)getProbe(x,y,z,a),stim,responseType,stateFlag,aftFlag,'UniformOutput', true);

nanFlag = isnan(probeBias + befBias + aftBias);
figure; subplot(1,2,1);
plotBar(befAcc(~nanFlag),probeAcc(~nanFlag),aftAcc(~nanFlag),'left');
subplot(1,2,2);
plotBar(abs(befBias(~nanFlag)),abs(probeBias(~nanFlag)),abs(aftBias(~nanFlag)),'right');


function plotBar(tempBefFlat,tempProbeFlat,tempAftFlat,tail)
    [hBef,pBef] = ttest(tempBefFlat,tempProbeFlat,'tail',tail);
    [hAft,pAft] = ttest(tempAftFlat,tempProbeFlat,'tail',tail);
    bar([nanmean(tempBefFlat) nanmean(tempProbeFlat) nanmean(tempAftFlat)] ,'EdgeColor','None','FaceColor',matlabColors(2,0.9)); hold on;
    for i = 1:length(tempBefFlat)
        f = plot([1 2 3],[tempBefFlat(i) tempProbeFlat(i) tempAftFlat(i)],'Color',[0.6 0.6 0.6],'Marker','.','MarkerSize',15,...
            'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',0.5);
    end
    legend(f,['pBef = ' num2str(pBef,'%.2e') newline 'pAft = ' num2str(pAft,'%.2e')],'Location','Best')
    xticks([1 2 3]); xticklabels({'Bef','Probe','Aft'}); xlim([0 4]);
end
    
end

function [acc,bias] = getProbe(stim,responseType,stateFlag,trialFlag)
    
    [~,acc_L,acc_R] = fn_getBias(stim,responseType,100);
    acc = acc_L/2 + acc_R/2; 
    accStart = 0.6; accEnd = max(acc)*0.85;
    startFlag = find(acc>accStart,1);
    endFlag = find(acc>accEnd,1)-1;

    reinfAccFlag = false(size(acc)); reinfAccFlag(startFlag:endFlag) = true;
    trialSelFlag = trialFlag & (stateFlag~=0) & reinfAccFlag;
    
    tempStim = stim(trialSelFlag); tempResp = responseType(trialSelFlag);
    if length(tempStim) > 20
        acc1 = nansum((tempStim==1) & (tempResp==1))/nansum(tempStim==1);
        acc2 = nansum((tempStim==2) & (tempResp==1))/nansum(tempStim==2);
        acc = acc1/2 + acc2/2; bias = acc1 - acc2;
    else; acc = nan; bias = nan;
    end
end

