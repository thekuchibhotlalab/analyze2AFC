%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);

%% Visualize zz068 day 5
mouse = 'zz068'; day = 5; plotBef = 20; plotAft = 10;
plotMouse = allMouse{strcmp(mouse,mice)};
probeSelIdx = find(plotMouse.behav.day==day & plotMouse.behav.probe==1);

tempTrialsIdx = probeSelIdx(1)-plotBef : probeSelIdx(1)+plotAft-1;
figure; subplot(2,1,1);
plotStuff(plotMouse,tempTrialsIdx); ylim([-1 4]);plot([plotBef plotBef] + 0.5,[-1 4],'Color',[0.8 0.8 0.8])

plotBef = 20; plotAft = 8;
tempTrialsIdx = probeSelIdx(11)-plotBef : probeSelIdx(11)+plotAft-1;
subplot(2,1,2);
plotStuff(plotMouse,tempTrialsIdx); ylim([-1 4]);plot([plotBef plotBef] + 0.5,[-1 4],'Color',[0.8 0.8 0.8])

%% Visualize zz068 day 3
mouse = 'zz068'; day = 3; plotBef = 20; plotAft = 8;
plotMouse = allMouse{strcmp(mouse,mice)};

probeSelIdx = find(plotMouse.behav.day==day & plotMouse.behav.probe==1);

tempTrialsIdx = probeSelIdx(1)-plotBef : probeSelIdx(1)+plotAft-1;
figure; subplot(2,1,1);
plotStuff(plotMouse,tempTrialsIdx); ylim([-1 4]);plot([plotBef plotBef] + 0.5,[-1 4],'Color',[0.8 0.8 0.8])

plotBef = 20; plotAft = 8;
tempTrialsIdx = probeSelIdx(9)-plotBef : probeSelIdx(9)+plotAft-1;
subplot(2,1,2);
plotStuff(plotMouse,tempTrialsIdx); ylim([-1 4]);plot([plotBef plotBef] + 0.5,[-1 4],'Color',[0.8 0.8 0.8])

%% Visualize zz107 day 4, wheel trajectories
mouse = 'zz107'; day = 4; plotBef = 20; plotAft = 9;
plotMouse = allMouse{strcmp(mouse,mice)};

probeSelIdx = find(plotMouse.behav.day==day & plotMouse.behav.probe==1);

tempTrialsIdx = probeSelIdx(1)-plotBef : probeSelIdx(1)+plotAft-1;
figure; %subplot(2,1,1);
plotStuff(plotMouse,tempTrialsIdx); ylim([-1 4]);plot([plotBef plotBef] + 0.5,[-1 4],'Color',[0.8 0.8 0.8])

title(['Reinf acc = 0.72, bias = -0.55' newline 'Probe acc = 1.0, bias = 0.0'])

%plotBef = 20; plotAft = 9;
%tempTrialsIdx = probeSelIdx(10)-plotBef : probeSelIdx(10)+plotAft-1;
%subplot(2,1,2);
%plotStuff(plotMouse,tempTrialsIdx); ylim([-1 4]);plot([plotBef plotBef] + 0.5,[-1 4],'Color',[0.8 0.8 0.8])
fn_plotWheelByTrial(mouseMega,mouse,postLearningTrial{exampleAnimal}+50);
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
    scatter(trialNum,action,25,colors,'filled')

    scatter(trialNum(tempResp == 2),stim(tempResp == 2),25,fn_wheelColorsPT('correct'))
    %scatter(trialNum(flagR),action,10,matlabColors(2),'filled')
    %legend('Stim 1','Stim 2','AutoUpdate','off')
    yticks([1 2]);yticklabels({'L','R'});
    xlabel('Trials')
end 
