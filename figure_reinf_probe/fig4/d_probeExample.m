%% LOAD ANIMALS
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);
reinf = mouseMega.loadReinf;

%% Visualize zz107 day 3, wheel trajectories, BLOCK 1
mouse = 'zz107'; day = 3; plotBef = 20; plotAft = 20; nProbe = 9; 
plotMouse = allMouse{strcmp(mouse,mice)};

probeSelIdx = find(plotMouse.behav.day==day & plotMouse.behav.probe==1);

tempTrialsIdx = probeSelIdx(1)-plotBef : probeSelIdx(1)+nProbe+plotAft-1;
figure; %subplot(2,1,1);
plotStuff(plotMouse,tempTrialsIdx);plot([0.6 2.4],[plotBef plotBef] + 0.5,'Color',[0.8 0.8 0.8])

title(['Reinf acc = 0.72, bias = -0.55' newline 'Probe acc = 1.0, bias = 0.0'])

fn_plotWheelByTrial(mouseMega,8,tempTrialsIdx,'nocolor','probeIdx',plotBef+1:plotBef+nProbe);

%% Visualize zz107 day 3, wheel trajectories, BLOCK 2
mouse = 'zz107'; day = 3; plotBef = 20; plotAft = 23; nProbe = 7; 
plotMouse = allMouse{strcmp(mouse,mice)};

probeSelIdx = find(plotMouse.behav.day==day & plotMouse.behav.probe==1);
tempTrialsIdxProbe = probeSelIdx(10):probeSelIdx(10)+nProbe-1;
tempTrialsIdx = probeSelIdx(10)-plotBef : probeSelIdx(10)+nProbe+plotAft-1;
tempTrialsIdxReinf = [probeSelIdx(10)-plotBef : probeSelIdx(10)-1 probeSelIdx(10)+nProbe: probeSelIdx(10)+nProbe+plotAft-1];
figure; %subplot(2,1,1);
plotStuff(plotMouse,tempTrialsIdx);plot([0.6 2.4],[plotBef plotBef] + 0.5,'Color',[0.8 0.8 0.8])

[bias, acc] = fn_getAccBias(plotMouse.behav.stimulus(tempTrialsIdxReinf), ...
    plotMouse.behav.stimulus(tempTrialsIdxReinf)==plotMouse.behav.action(tempTrialsIdxReinf));
[biasP, accP] = fn_getAccBias(plotMouse.behav.stimulus(tempTrialsIdxProbe), ...
    plotMouse.behav.stimulus(tempTrialsIdxProbe)==plotMouse.behav.action(tempTrialsIdxProbe));
title(['Reinf acc = ' num2str(acc) ', bias =' num2str(abs(bias)) newline 'Probe acc = ' num2str(accP) ', bias = ' num2str(abs(biasP))])

fn_plotWheelByTrial(mouseMega,8,tempTrialsIdx,'nocolor','probeIdx',plotBef+1:plotBef+nProbe);



%% Visualize zz107 day 4, wheel trajectories, BLOCK 1
mouse = 'zz107'; day = 4; plotBef = 20; plotAft = 22; nProbe = 8; 
plotMouse = allMouse{strcmp(mouse,mice)};

probeSelIdx = find(plotMouse.behav.day==day & plotMouse.behav.probe==1);

tempTrialsIdx = probeSelIdx(1)-plotBef : probeSelIdx(1)+nProbe+plotAft-1;
figure; %subplot(2,1,1);
plotStuff(plotMouse,tempTrialsIdx);plot([0.6 2.4],[plotBef plotBef] + 0.5,'Color',[0.8 0.8 0.8])

title(['Reinf acc = 0.72, bias = -0.55' newline 'Probe acc = 1.0, bias = 0.0'])

fn_plotWheelByTrial(mouseMega,8,tempTrialsIdx,'nocolor','probeIdx',plotBef+1:plotBef+nProbe);

%% Visualize zz107 day 4, wheel trajectories, BLOCK 2
mouse = 'zz107'; day = 4; plotBef = 20; plotAft = 20; nProbe = 10; 
plotMouse = allMouse{strcmp(mouse,mice)};

probeSelIdx = find(plotMouse.behav.day==day & plotMouse.behav.probe==1);

tempTrialsIdx = probeSelIdx(9)-plotBef : probeSelIdx(9)+nProbe+plotAft-1;
figure; %subplot(2,1,1);
plotStuff(plotMouse,tempTrialsIdx);plot([0.6 2.4],[plotBef plotBef] + 0.5,'Color',[0.8 0.8 0.8])

title(['Reinf acc = 0.72, bias = -0.55' newline 'Probe acc = 1.0, bias = 0.0'])

fn_plotWheelByTrial(mouseMega,8,tempTrialsIdx,'nocolor','probeIdx',plotBef+1:plotBef+nProbe);

%% functions
function plotStuff(plotMouse,tempTrialsIdx)
    stim = plotMouse.behav.stimulus(tempTrialsIdx);% flagR = plotMouse.behav.stimulus(tempTrialsIdx) == 2;
    action = plotMouse.behav.action(tempTrialsIdx);
    trialNum = 1:length(tempTrialsIdx);
    colors = zeros(length(trialNum),3);
    tempResp = plotMouse.behav.responseType(tempTrialsIdx);
    colors(tempResp == 1,:) = repmat(fn_wheelColorsPT('correct'),[sum(tempResp == 1) 1]);
    colors(tempResp == 2,:) = repmat(fn_wheelColorsPT('incorrect'),[sum(tempResp == 2) 1]);
    hold on;
    scatter(action,trialNum,25,colors,'filled')

    scatter(stim(tempResp == 2),trialNum(tempResp == 2),25,fn_wheelColorsPT('correct'))
    %scatter(trialNum(flagR),action,10,matlabColors(2),'filled')
    %legend('Stim 1','Stim 2','AutoUpdate','off')
    xticks([1 2]);xticklabels({'L','R'});
    ylabel('Trials'); 
end 