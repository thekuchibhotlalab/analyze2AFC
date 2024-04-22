%% GET OBJECT FOR ALL THE PUREOTNE ANIMALS 
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone';

allMouse = fn_getObjPT(mice,task);



%%  GET OBJECT FOR OTHER COHORTS
clear; 
mice ={'zz101','zz107'};
task = {'PT_pip_task2','FM_Dir_task2'};

allMouse = getObj(mice,task);

%%  GET OBJECT FOR FM Dir Dur COHORTS
clear; 
mice ={'zz101','zz102','zz105'};
task = 'FM_Dir_Dur_task1';

allMouse = getObj(mice,task);  

%%  GET OBJECT FOR FM Dir Dur COHORTS
clear; 
mice ={'zz121','zz122','zz123'};
task = 'FM_Dir_task1';

allMouse = getObj(mice,task);  
%%  GET OBJECT FOR PT pip COHORTS
clear; 
mice ={'zz097','zz098'};
task = {'PT_pip'};

allMouse = getObj(mice,task);   

%%  GET OBJECT FOR PT pip COHORTS
clear; 
mice ={'zz097','zz098','msb01','msb02','msb03','msb04'};
task = {'PT_pip'};

allMouse = getObj(mice,task); 

plotPerfBehavLR (allMouse{1});
plotPerfBehavLR (allMouse{2});
plotPerfBehavLR (allMouse{3});
plotPerfBehavLR (allMouse{4});
%%
cellfun(@plotPerfBehavIndiv,allMouse);
%%
cellfun(@plotBlockTransition,allMouse);

%% MouseMEGA plot
mouseMega = wheel2AFCmega(allMouse);

%%
plotBlockTransition(mouseMega);
%%
plotBlockTransitionFreq(mouseMega)
%% 
plotReinfProbe(mouseMega);
%%
plotBiasTransitions(mouseMega);
%%
plotReactionTime(mouseMega);

%% MAKE PROBE IN BIAS VS. NOT BIAS BLOCK PLOT
mice = {'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz097','zz098','zz107','zz109'};
rootPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\';
loadPath = [rootPath filesep 'trialData\'];
allMouse = {};
for i = 1:length(mice)
    load([loadPath mice{i} '.mat']);
    tempObj = wheel2AFC(trialData,'opsPath','C:\Users\zzhu34\Documents\tempdata\octoData\opsFiles','mouse',mice{i} );
    tempObj = tempObj.removeMiss(); tempObj = tempObj.calculateProbe();
    tempObj = tempObj.getAcc();
    tempObj = tempObj.getBiasBlock(0.3);
    allMouse{i} = tempObj;
end
clear rootPath loadPath tempObj i trialData
mouseMega = wheel2AFCmega(allMouse);
plotReinfProbeHighBias(mouseMega);

%% ---------------------PLOT nBias Blocks----------------------
blockLen = []; 
for i = 1:length(allMouse)
    blockLen = cat(1,blockLen,[allMouse{i}.biasBlock.blockL.len;allMouse{i}.biasBlock.blockR.len]);
    nBlockL(i) = length(allMouse{i}.biasBlock.blockL.len);
    nBlockR(i) = length(allMouse{i}.biasBlock.blockR.len);
end
temp = cellfun(@(x)num2str(x,'%02d'),num2cell(nBlockL),'UniformOutput',false);
disp(['blockL: ' strjoin(temp,' ')]); 
temp = cellfun(@(x)num2str(x,'%02d'),num2cell(nBlockR),'UniformOutput',false);
disp(['blockR: ' strjoin(temp,' ')]); 

figure; histogram(blockLen,0:2:40); ylabel('block length'); xlabel('frequency')

%% ---------------------PLOT ACCURACY AND BIAS TOGETHER----------------------
fn_figureSmartDim('hSize',0.3,'widthHeightRatio',1) ; hold on; 
avgAccSmooth = smoothdata(nanmean(accuracy,2),'movmean',200);
avgBiasSmooth = smoothdata(nanmean(bias,2),'movmean',200);
for i = 1:size(accuracy,2)
    accSmooth = smoothdata(accuracy(:,i),'movmean',400); biasSmooth = smoothdata(bias(:,i),'movmean',400);
    plot(accSmooth,biasSmooth,'Color',matlabColors(i,0.3));
    plot(accSmooth(1),biasSmooth(1),'.','MarkerSize',15,'Color',matlabColors(i))
    plot(accSmooth(end),biasSmooth(end),'*','MarkerSize',8,'Color',matlabColors(i))
end
plot(avgAccSmooth,avgBiasSmooth,'Color',[0 0 0],'LineWidth',2);
xlim([0.3 0.95]); xticks([0.3 0.5 0.7 0.9])
xlabel('Accuracy'); ylabel('Action Bias')

%% ------------------------------REACTION TIME-------------------------------
reactionTime = mouseMega.reactionTime;
reactionTime(reactionTime > 2.6) = nan; reactionTime(reactionTime < 0) = nan;
reactionTime = smoothdata(reactionTime,1,'movmean',50,'includenan');

fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
for j = 1:size(reactionTime,2)
    plot(reactionTime(:,j),'Color', matlabColors(1,0.2), 'LineWidth',1);
end
plot(nanmean(reactionTime,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
ylim([0.3 1.5]); xlim([1 trialLim]);
ylabel('nMovements'); xlabel('Trials');

%% ------------------------------ITI MOVEMENTS-------------------------------
wheelPreSound = mouseMega.wheel.wheelPreSound;
wheelPreSound = smoothdata(wheelPreSound,1,'movmean',50,'includenan');
%wheelPreSound = wheelPreSound./repmat(nanmean(wheelPreSound(1:500,:),1),trialLim,1); %normalize to the first day of each animal
fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); hold on;
for j = 1:size(accuracy,2)
    plot(wheelPreSound(:,j),'Color', matlabColors(1,0.2), 'LineWidth',1);
end
plot(nanmean(wheelPreSound,2),'Color', matlabColors(1,0.9), 'LineWidth',2);
xlim([1 trialLim]);
ylabel('Reaction Time (s)'); xlabel('Trials');

%% ------------------------------FUNCTIONS-------------------------------
function allMouse = getObj(mice,task)
    
rootPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\';
loadPath = [rootPath filesep 'trialData\'];

allMouse = {};
for i = 1:length(mice)
    load([loadPath mice{i} '.mat']);

    

    ops.datapath = loadPath;
    ops.behavType = 'wheel2AFC';

    ops.preprocessedInput = false;
    ops.psytrackPath = [ rootPath '\psyTrackData\psyTrackFit\'];
    ops.correctWheelTrial = true;
    ops.learningCurveBin = 300;
    ops.selectProtocol = task;
    tempObj = wheel2AFC(trialData,'ops',ops,'mouse',mice{i} );
    
    
    tempObj = tempObj.removeMiss(); 
    %tempObj = tempObj.removeProbe();  
    tempObj = tempObj.getAcc();
    tempObj = tempObj.calculateProbe();
    tempObj = tempObj.getBiasBlock(0.2);
    allMouse{i} = tempObj;
end
    
    
end
