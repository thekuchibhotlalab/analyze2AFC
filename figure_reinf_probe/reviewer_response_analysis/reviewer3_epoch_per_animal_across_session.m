%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT_bin40(mice);
learningCurveBin = 40;
mouseMega = wheel2AFCmega(allMouse);

%%

plotBlockTransition(mouseMega)

%% accuray of each day
figure;
for i = 1:13; temp = []; subplot(3,5,i);
    for j = 1:max(allMouse{i}.behav.day)
        tempFlag = allMouse{i}.behav.day==j; 
        [bias, acc] = fn_getAccBias(allMouse{i}.behav.stimulus(tempFlag),...
            allMouse{i}.behav.stimulus(tempFlag)== allMouse{i}.behav.action(tempFlag));
    temp(j) = acc;
    end
    plot(temp);
end