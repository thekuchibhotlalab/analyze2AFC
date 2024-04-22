%% LOAD ANIMALS
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
allMouse = fn_getObjPT_bin40(mice);
mouseMega = wheel2AFCmega(allMouse);
reinf = mouseMega.loadReinf;

nMouse = size(reinf.acc,2);
plotLim = 3000;


%%
maxAcc = cellfun(@(x)(max(x.behav.acc)),allMouse,'UniformOutput',true);
maxAccModel = cellfun(@(x)(max(x.behav.modelAcc)),allMouse,'UniformOutput',true);
maxAccTrial = []; acc60Trial = [];acc70Trial = [];acc80Trial = []; stimCorr = []; rewardRate = [];
for i = 1:length(allMouse)
    [~,idx]= max(allMouse{i}.behav.acc); maxAccTrial(i)=idx;
    [idx]= find(allMouse{i}.behav.modelAcc>0.8,1,'first'); acc80Trial(i)=idx;
    [idx]= find(allMouse{i}.behav.modelAcc>0.7,1,'first'); acc70Trial(i)=idx;
    [idx]= find(allMouse{i}.behav.modelAcc>0.5,1,'first'); acc60Trial(i)=idx;
    meanBiasMid(i)= nanmean(abs(allMouse{i}.behav.bias(acc60Trial(i):acc80Trial(i))));
    stimCorr(i)= abs(mean(allMouse{i}.behav.stimulus(acc60Trial(i):acc80Trial(i))==1)-...
        mean(allMouse{i}.behav.stimulus(acc60Trial(i):acc80Trial(i))==2));
    rewardRate(i)= mean(allMouse{i}.behav.reward(acc60Trial(i):acc80Trial(i)));
end 
 
meanBias = cellfun(@(x)(nanmean(abs(x.behav.bias))),allMouse,'UniformOutput',true);

%figure; scatter(meanBias,maxAccTrial);
%figure; fn_plotScatterCorr(meanBias,acc70Trial);
%figure; scatter(meanBias,maxAcc);
tempColors = zeros(length(meanBiasMid),3);
for i = 1:7; tempColors(i,:) = matlabColors(3); end 
for i = 8:13; tempColors(i,:) = matlabColors(4); end 
figure; [lm,sc] = fn_plotScatterCorr(meanBiasMid,acc80Trial-acc60Trial,'scatterSize',30,'scatterColors', {tempColors}, 'lineColor',{[0.2 0.2 0.2]});
xlim([0 0.8]); ylim([0 4000]);
xlabel('avg. bias magnitude'); ylabel(['Trials to take from 0.55 accuracy' newline 'to 0.75 accuracy (model predicted)'])

learningSpeed = acc80Trial-acc60Trial;
figure; 
[p] = fn_plotComparison({learningSpeed(1:7), learningSpeed(8:13)});
xticks([1 2]); xticklabels({'2sec timeout','5 sec timeout'})
ylabel(['Trials to take from 0.6 accuracy' newline 'to 0.8 accuracy (model predicted)'])

%figure; [lm,sc] = fn_plotScatterCorr(meanBiasMid,stimCorr,'scatterSize',30,'scatterColors', {tempColors}, 'lineColor',{[0.2 0.2 0.2]});
%xlabel('avg. bias magnitude'); ylabel('stim correction -- abs(L freq - R freq)')
%xlim([0 0.6]); ylim([0 0.4])

%figure; [lm,sc] = fn_plotScatterCorr(meanBiasMid,rewardRate,'scatterSize',30,'scatterColors', {tempColors}, 'lineColor',{[0.2 0.2 0.2]});
%xlabel('avg. bias magnitude'); ylabel('stim correction -- abs(L freq - R freq)')
%xlim([0 0.6]); ylim([0.4 0.8])


%%