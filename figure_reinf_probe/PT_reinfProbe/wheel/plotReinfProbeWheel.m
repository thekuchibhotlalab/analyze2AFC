%%
clear; 
mice ={'zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone';

rootPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\';
loadPath = [rootPath filesep 'trialData\'];

allMouse = {};
for i = 1:length(mice)
    load([loadPath mice{i} '.mat']);
    if iscell(task)
        tempObj = wheel2AFC(trialData,'opsPath',[rootPath filesep 'opsFiles' filesep task{i} ],'mouse',mice{i} );
    else
        tempObj = wheel2AFC(trialData,'opsPath',[rootPath filesep 'opsFiles' filesep task ],'mouse',mice{i} );
    end
    tempObj = tempObj.removeMiss(); 
    tempObj = tempObj.getAcc();
    tempObj = tempObj.calculateProbe();
    allMouse{i} = tempObj;
end

%%
f3 = figure;
for i = 1:length(mice)
    f1 = figure; f2 = figure; 
    nDay = unique(allMouse{i}.behav.day);
    rtReinfBef = [];rtProbe = [];
    for j = 1:length(nDay)
        dayFlag = allMouse{i}.behav.day==nDay(j);
        probeFlag = allMouse{i}.behav.goodProbe==1;
        reinfBefFlag = allMouse{i}.behav.reinfBef==1;

        reinfWheel = allMouse{i}.behavVar.wheelPos_aligned(~probeFlag & dayFlag,:);
        reinfBefWheel = allMouse{i}.behavVar.wheelPos_aligned(reinfBefFlag & dayFlag,:);
        probeWheel = allMouse{i}.behavVar.wheelPos_aligned(probeFlag & dayFlag,:);
        [a,b] = fn_sqrtInt(length(nDay));
        figure(f1); subplot(a,b,j); hold on; 
        plot(nanmean(abs(reinfWheel),1)); plot(nanmean(abs(probeWheel),1)); plot(nanmean(abs(reinfBefWheel),1)); 
        figure(f2); subplot(a,b,j); hold on; 
        for k = 1:size(probeWheel,1); plot(abs(probeWheel(k,:)),'Color',matlabColors(2)); end 
        for k = 1:size(reinfBefWheel,1); plot(abs(reinfBefWheel(k,:)),'Color',matlabColors(3)); end 
        
        rtReinfBef(j) = mean(allMouse{i}.behav.reactionTime(reinfBefFlag & dayFlag,:));
        rtProbe(j) = mean(allMouse{i}.behav.reactionTime(probeFlag & dayFlag,:));
    end
    figure(f3); [a,b] = fn_sqrtInt(length(mice));
    subplot(a,b,i); hold on; plot(rtReinfBef); hold on; plot(rtProbe);
end

%%
selDay  = [1 4 8];
mouse = 6; figure; 
for i = 1:length(selDay)
    nDay = unique(allMouse{mouse}.behav.day);
    dayFlag = allMouse{mouse}.behav.day==nDay(selDay(i));
    probeFlag = allMouse{mouse}.behav.context==3;
    reinfWheel = allMouse{mouse}.behavVar.wheelPos_aligned(~probeFlag & dayFlag,:);

    subplot(2,length(selDay),i); plot(reinfWheel(1:30,:)'); xlim([1000 4501])
    subplot(2,length(selDay),length(selDay)+i); plot(reinfWheel(end-30:end,:)'); xlim([1000 4501])
end