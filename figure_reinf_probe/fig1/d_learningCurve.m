%% LOAD ANIMALS
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);

%% PLOT 1 -- Learning curve of accuracy, bias, action rate and RT
reinf = mouseMega.loadReinf;
reinf.reactionTime(reinf.reactionTime>=2.5 | reinf.reactionTime < 0) = nan;
reinf.reactionTime = smoothdata(reinf.reactionTime,1,'movmean',30);

nMouse = size(reinf.acc,2);
plotLim = 3000;

fn_figureSmartDim('hSize',0.6,'widthHeightRatio',0.4); 
subplot(3,1,1);hold on;
plot([1 plotLim],[0.5 0.5],'Color',[0.6 0.6 0.6],'LineWidth',2);xlim([1 plotLim]);
plotIndividualReinf(reinf.acc(1:plotLim,:));ylabel('Accuracy'); ylim([0.3 1.0]); yticks([0.5 0.75 1])


subplot(3,1,2);hold on;
plotIndividualReinf(reinf.bias(1:plotLim,:));
xlim([1 plotLim]);ylabel('Bias'); ylim([0 1.0])
xlabel('Trials');

subplot(3,1,3);hold on;
plotIndividualReinf(smoothdata(reinf.reactionTime(1:plotLim,:),1,'movmean',100,'includenan'));ylabel('reactionTime'); ylim([0 1.8])
xlim([1 plotLim]); xlabel('Trials'); ylim([0 1.5]); yticks(0:0.5:1.5)

%% PLOT 2 -- Naive vs. Expert acc and RT
% select 100 trials for acc and RT comparison
reinf.reactionTime(reinf.reactionTime>=2.5 | reinf.reactionTime < 0) = nan;
reinf.reactionTime = smoothdata(reinf.reactionTime,1,'movmean',30);

tempBin = 300; start = 1;
preLearningTrial = start:start+tempBin-1;
[~,accPre,dpPre] = cellfun(@(x)(fn_getAccBias(x.behav.stimulus(preLearningTrial),x.behav.responseType(preLearningTrial)==1)),...
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
% get mean acc and RT
[~,accPost,dpPost] = cellfun(@(x,y)(fn_getAccBias(x.behav.stimulus(y),x.behav.responseType(y)==1)),...
    mouseMega.mouseCell,postLearningTrial);
rtPre = cellfun(@(x)(mean(x.behav.reactionTime(preLearningTrial))),...
    mouseMega.mouseCell);
rtPost = cellfun(@(x,y)(mean(x.behav.reactionTime(y))),...
    mouseMega.mouseCell,postLearningTrial);
% plot results
fn_figureSmartDim('hSize',0.7,'widthHeightRatio',0.6); 
subplot(3,1,1);
[p] = fn_plotComparison({accPre,accPost},'paired',true,'dotType','side','compType','errorbarWithDot','scatterArgIn',...
    {5,'MarkerEdgeColor','none','MarkerFaceColor','none'},'errorbarArgIn', {'Color',[0.2 0.2 0.2],'LineWidth',1.5,'LineStyle','none'}); 
xlim([0.6 2.4]);  ylim([0.25 1]); yticks([0.25 0.5 0.75 1]); title(['acc, p=' num2str(p)])

subplot(3,1,2);
[p] = fn_plotComparison({rtPre,rtPost},'paired',true,'dotType','side','compType','errorbarWithDot','scatterArgIn',...
    {5,'MarkerEdgeColor','none','MarkerFaceColor','none'},'errorbarArgIn', {'Color',[0.2 0.2 0.2],'LineWidth',1.5,'LineStyle','none'}); 
xlim([0.6 2.4]);  ylim([0 1.5]); yticks(0:0.5:1.5); title(['RT, p=' num2str(p)])

% redefine the 
tempBin = 300;
preLearningTrial = {};
[~,startIdx] = max(abs(reinf.bias(1:1000,:)),[],1);
for i = 1:mouseMega.nMouse
    if startIdx(i) - tempBin/2 <= 0
        preLearningTrial{i} = 1:tempBin;
    else
        preLearningTrial{i} = startIdx(i) - tempBin/2 + 1 : startIdx(i) + tempBin/2;
    end  
end
biasPre = cellfun(@(x,y)(nanmean(abs(x.behav.bias(y)))),...
    mouseMega.mouseCell,preLearningTrial);
biasPost = cellfun(@(x,y)(nanmean(abs(x.behav.bias(y)))),...
    mouseMega.mouseCell,postLearningTrial);

subplot(3,1,3);
[p] = fn_plotComparison({abs(biasPre),abs(biasPost)},'paired',true,'dotType','side','compType','errorbarWithDot','scatterArgIn',...
    {5,'MarkerEdgeColor','none','MarkerFaceColor','none'},'errorbarArgIn', {'Color',[0.2 0.2 0.2],'LineWidth',1.5,'LineStyle','none'}); 
xlim([0.6 2.4]);  ylim([0 1]); yticks([0 0.5 1]); title(['bias, p=' num2str(p)])

%% FUNCTIONS
function plotIndividualReinf(mat)
    % PLOT REINF
    hold on;
    plot(mat,'Color', fn_wheelColorsPT('Reinf',0.3), 'LineWidth',1);
    plot(nanmean(mat,2),'Color', fn_wheelColorsPT('Reinf'), 'LineWidth',3);
end

