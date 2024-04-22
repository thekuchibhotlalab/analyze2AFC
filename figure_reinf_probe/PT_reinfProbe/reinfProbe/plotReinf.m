%%
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};

allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);

%%
clear; 
mice ={'zz121','zz122','zz123'};
task = 'FM_Dir_task1';
allMouse = fn_getObjPT(mice,task);
mouseMega = wheel2AFCmega(allMouse);
%%
reinf = mouseMega.loadReinf;
probeBin = mouseMega.loadProbeBin;
reinf.reactionTime(reinf.reactionTime>=2.5 | reinf.reactionTime < 0) = nan;
reinf.reactionTime = smoothdata(reinf.reactionTime,1,'movmean',30);
%reinf.reactionTime(reinf.reactionTime>=1.8 | reinf.reactionTime < 0) = nan;

nMouse = size(reinf.acc,2);
plotLim = 3000;

fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8); 
subplot(2,2,1);hold on;
plot([1 plotLim],[0.5 0.5],'Color',[0.6 0.6 0.6],'LineWidth',2);xlim([1 plotLim]);
plotIndividualReinf(reinf.acc(1:plotLim,:));ylabel('Accuracy'); ylim([0.3 1.0]); yticks([0.5 0.75 1])

subplot(2,2,3);hold on;
plotIndividualReinf(reinf.bias(1:plotLim,:));
xlim([1 plotLim]);ylabel('Bias'); ylim([0 0.8])
xlabel('Trials');
 
subplot(2,2,2);hold on;
plotIndividualReinf(reinf.actionRate(1:plotLim,:));
xlim([1 plotLim]); ylabel('actionRate'); ylim([0.2 1.0])

subplot(2,2,4);hold on;
plotIndividualReinf(smoothdata(reinf.reactionTime(1:plotLim,:),1,'movmean',100,'includenan'));ylabel('reactionTime'); ylim([0 1.8])
xlim([1 plotLim]); xlabel('Trials'); ylim([0 1.5]); yticks(0:0.5:1.5)
%% plot naive vs. expert animal
tempBin = 100; start = 100;
preLearningTrial = start:start+tempBin-1;
[~,~,dpPre] = cellfun(@(x)(fn_getAccBias(x.behav.stimulus(preLearningTrial),x.behav.responseType(preLearningTrial)==1)),...
    mouseMega.mouseCell);

rtPre = cellfun(@(x)(mean(x.behav.reactionTime(preLearningTrial))),...
    mouseMega.mouseCell);

[~,endIdx] = max(reinf.acc,[],1);
postLearningTrial = {};
for i = 1:mouseMega.nMouse
    if endIdx(i) + tempBin/2 > mouseMega.nTrials(i)
        postLearningTrial{i} = (mouseMega.nTrials(i)-tempBin+1): mouseMega.nTrials(i);
    else
        postLearningTrial{i} = endIdx(i) - tempBin/2 + 1 : endIdx(i) + tempBin/2;
    end  
end

[~,~,dpPost] = cellfun(@(x,y)(fn_getAccBias(x.behav.stimulus(y),x.behav.responseType(y)==1)),...
    mouseMega.mouseCell,postLearningTrial);
rtPost = cellfun(@(x,y)(mean(x.behav.reactionTime(y))),...
    mouseMega.mouseCell,postLearningTrial);


figure;
[p] = fn_plotBarPaired({rtPre,rtPost});


%%

[wheelDownSample, RT] = loadWheelClusterData(mouseMega, 8:13);
exampleAnimal = 10;
% post learning
fn_plotWheelByStim(mouseMega,exampleAnimal,postLearningTrial{exampleAnimal}+50)
% pre learning
fn_plotWheelByStim(mouseMega,exampleAnimal,preLearningTrial-50)
%figure;
%fn_plotWheelByTrial(mouseMega,exampleAnimal,postLearningTrial{exampleAnimal}+50);
%% functions


function plotIndividualReinf(mat)
    % PLOT REINF
    plot(mat,'Color', fn_wheelColorsPT('Reinf',0.3), 'LineWidth',1);
    plot(nanmean(mat,2),'Color', fn_wheelColorsPT('Reinf'), 'LineWidth',3);
end

function [wheelDownSample, RT] = loadWheelClusterData(mouseMega, idx)
    if isempty(idx); idx = 1:mouseMega.nMouse; end
    nTrialsum = [0 cumsum(mouseMega.nTrials)];
    [wheelDownSample,badFlag] = getDownSample(mouseMega,idx);
    RT = mouseMega.getProp('behav','field','reactionTime','matFlag',false,'idx',idx);
end

function [wheelDownSample,badFlag] = getDownSample(mouseMega,idx)
if isempty(idx); idx = 1:mouseMega.nMouse; end
wheelDownSample = cell(1,length(idx));
for i = 1:length(idx)
    [wheelPos, badFlag] = removeBadTrials(mouseMega.mouseCell{idx(i)}.behavVar.wheelPos_aligned);
    
    % normalize the max position of the wheel to 1; downsample every 50ms 
    downSampleRate = 50; tempTimeDownSampleIdx = 2001:downSampleRate:4501;
    tempMax = wheelPos(:,tempTimeDownSampleIdx(1):tempTimeDownSampleIdx(end)); tempMax = max(tempMax(:));
    wheelPos = wheelPos /tempMax;
    tempWheelDownSample = wheelPos(:,tempTimeDownSampleIdx);
    wheelDownSample{i} = tempWheelDownSample;
end 
end

function [tempWheel, badFlag] = removeBadTrials(tempWheel)
    totalTrial = size(tempWheel,1);
    badFlag = tempWheel(:,end)==0;
    tempWheel(badFlag,:) = nan;
    disp(['Bad Trials = ' int2str(sum(badFlag)) ' out of ' int2str(totalTrial)])
end