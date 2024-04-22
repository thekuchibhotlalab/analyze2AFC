%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT_bin30(mice);
for i = 1:length(mice); allMouse{i} = allMouse{i}.getPsyTrack(); end

mouseMega = wheel2AFCmega(allMouse);

%% find the bin that best correspond
learningCurveBin = 20:10:500;
tempCorr = zeros(length(mice),length(learningCurveBin));
for i = 1:length(mice)
    for j = 1:length(learningCurveBin)
        [tempBias] = fn_getAccBiasSmooth(allMouse{i}.behav.stimulus,allMouse{i}.behav.responseType,learningCurveBin(j));
        tempCorr(i,j) = corr(tempBias, allMouse{i}.behav.modelBias,'rows','complete');
    end 
end

%%
[binVal,binIdx] = max(tempCorr,[],2);
exampeleAni = 11;
figure; subplot(2,1,1);
[tempBias] = fn_getAccBiasSmooth(allMouse{exampeleAni}.behav.stimulus,allMouse{exampeleAni}.behav.responseType,30);
plot(tempBias); hold on; plot(allMouse{exampeleAni}.behav.modelBias);
xlim([0 2500]);xticks(0:500:2500)
subplot(2,1,2);
[tempBias] = fn_getAccBiasSmooth(allMouse{exampeleAni}.behav.stimulus,allMouse{exampeleAni}.behav.responseType,learningCurveBin(binIdx(exampeleAni)));
plot(tempBias); hold on; plot(allMouse{exampeleAni}.behav.modelBias);
xlim([0 2500]);xticks(0:500:2500)

tempBest = [];
for i = 1:length(binIdx); tempBest(i) = learningCurveBin(binIdx(i)); end 
figure; subplot(2,1,1); cdfplot(tempBest);
subplot(2,1,2); cdfplot(binVal); xlim([0.86 1])



