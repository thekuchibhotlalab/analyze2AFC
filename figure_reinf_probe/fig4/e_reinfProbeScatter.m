%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);

probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;

%% PLOT 1.0 -- Subsample within DAY
probeDay = mouseMega.objFun('computeProbeAlignByTrial',{'day',[reinfThre probeThre],probeTrialBin});

%% PLOT 1.1 -- reinf vs. probe comparison
probeRand = mouseMega.objFun('computeProbeAlignByTrialRandomize',{}); 

c1 = fn_cell2mat(probeRand{1},2);

lm = fitlm(c1(1,:),c1(2,:));%figure; scatter(c1(1,:),c1(2,:),20,'filled')
c2 = fn_cell2mat(probeRand{2},2);
%lm2 = fitlm(c2(1,:),c2(2,:));%figure; scatter((c2(1,:)),c2(2,:),20,'filled')
plotSubsampleAnimal(probeDay,'randomComp',{c1,c2});
%% PLOT 1.2 -- divide probe blocks into high vs. low violations
probeDay = mouseMega.objFun('computeProbeAlignByTrial',{'befAft',[reinfThre probeThre],probeTrialBin});

%% FUNCTIONS
function plotSubsampleAnimal(probeCount,varargin)
p = inputParser;p.KeepUnmatched = true;
p.addParameter('randomComp', {}); p.parse(varargin{:});

% plot by animal
probeAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,1)'))),probeCount{2},'UniformOutput',false)); 
probeBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.probe(:,2)'))),probeCount{2},'UniformOutput',false)); 
reinfAcc = fn_cell2matFillNan(cellfun(@(x)((abs(x.reinf(:,1)'))),probeCount{2},'UniformOutput',false)); 
reinfBias = fn_cell2matFillNan(cellfun(@(x)((abs(x.reinf(:,2)'))),probeCount{2},'UniformOutput',false)); 

colors = parula; tempIdx = round(linspace(1, size(colors,1),size(reinfBias,1)+1));
colors = colors(tempIdx(1:end-1),:);
tempColor = []; tempColor(:,1,:) = colors; tempColor = repmat(tempColor,1,size(probeAcc,2),1);
tempColor = reshape(tempColor,size(tempColor,1)* size(tempColor,2),[]);

% PLOT FOR BIAS VS ACC CORRELATION
fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.9); hold on;
if ~isempty(p.Results.randomComp)
    temp = p.Results.randomComp{1}; temp2 = p.Results.randomComp{2}; temp(1,:) = temp2(1,:);
    lm2 = fitlm(temp(1,:),temp(2,:)); [r1,p1] = corrcoef(temp','rows','complete');
    tempX = [-0 0.85]; tempY = tempX * lm2.Coefficients{2,1} + lm2.Coefficients{1,1};
    plot(tempX,tempY,'Color',[0.6 0.6 0.6],'LineWidth',0.5);
    dispFlag = randi(size(temp,2),1,500); scatter(temp(1,dispFlag),temp(2,dispFlag),15,[0.9 0.9 0.9],'filled');

    tempAxis = 0:0.01:0.85; [Ypred,YCI] = predict(lm2, tempAxis');
    fill([tempAxis fliplr(tempAxis)], [YCI(:,1);flipud(YCI(:,2))]',[0.6 0.6 0.6],'LineStyle','None','faceAlpha',0.2);
    %plot_gaussian_ellipsoid([nanmean(temp(1,:)) nanmean(temp(2,:))], cov(temp','omitrows'),1,{'LineWidth',2,'Color',[0.8 0.8 0.8]});
    
end

tempB = reinfBias(:)-probeBias(:); tempB = abs(reinfBias(:)); 
tempA = probeAcc(:)-reinfAcc(:); 
lm = fitlm(tempB,tempA);tempX = [0 0.85]; tempY = tempX * lm.Coefficients{2,1} + lm.Coefficients{1,1};
if ~isempty(p.Results.randomComp)
    %plot([0 0],[-0.4 0.4],'LineWidth',0.5,'Color',[0.8 0.8 0.8]); 
    plot([0 1],[0 0],'LineWidth',0.5,'Color',[0.8 0.8 0.8]);
    [r2,p2] = corrcoef(tempB,tempA,'rows','complete');
    %scatter(reinfBias(:)-probeBias(:),probeAcc(:)-reinfAcc(:),30,matlabColors(1),'filled'); hold on;
    scatter(tempB,tempA,30,fn_wheelColorsPT('probe'),'filled'); hold on;
    plot(tempX,tempY,'Color',fn_wheelColorsPT('probe'),'LineWidth',0.5);
    tempAxis = 0:0.01:0.85; [Ypred,YCI] = predict(lm, tempAxis');
    fill([tempAxis fliplr(tempAxis)], [YCI(:,1);flipud(YCI(:,2))]',fn_wheelColorsPT('probe'),'LineStyle','None','faceAlpha',0.2);
    title(['corr = ' num2str(r2(1,2)) ', p = ' num2str(p2(1,2)) ', slope = ' num2str(lm.Coefficients.Estimate(2)) newline ...
        'corr = ' num2str(r1(1,2))  ', p = ' num2str(p1(1,2)) ', slope = ' num2str(lm2.Coefficients.Estimate(2))])
    %plot_gaussian_ellipsoid([nanmean(tempB) nanmean(tempA)], cov(tempB,tempA,'omitrows'),1,{'LineWidth',2,'Color',matlabColors(1,0.4)});

else
    scatter(reinfBias(:)-probeBias(:),probeAcc(:)-reinfAcc(:),30,tempColor,'filled'); hold on;
    plot(tempX,tempY,'Color',[0.6 0.6 0.6],'LineWidth',0.5);
    title(['corr = ' num2str(corr(tempB,tempA,'rows','complete')) ', p = ' num2str(lm.Coefficients.pValue(2))])
end
xticks(0:0.25:1); 
xlim([0 1]); ylim([-0.5 0.5]);yticks(-0.5:0.25:0.5);



% PLOT FOR BIAS VS BIAS REDUCTION CORRELATION
fn_figureSmartDim('hSize',0.5,'widthHeightRatio',0.9); hold on;
if ~isempty(p.Results.randomComp)
    temp = p.Results.randomComp{2}; 
    lm2 = fitlm(temp(1,:),temp(2,:));
    tempX = [0 0.85]; tempY = tempX * lm2.Coefficients{2,1} + lm2.Coefficients{1,1};
    [r1,p1] = corrcoef(temp','rows','complete');
    plot(tempX,tempY,'Color',[0.6 0.6 0.6],'LineWidth',0.5);
    dispFlag = randi(size(temp,2),1,500); scatter(temp(1,dispFlag),temp(2,dispFlag),15,[0.9 0.9 0.9],'filled');
    
    tempAxis = 0:0.01:0.85; [Ypred,YCI] = predict(lm2, tempAxis');
    fill([tempAxis fliplr(tempAxis)], [YCI(:,1);flipud(YCI(:,2))]',[0.6 0.6 0.6],'LineStyle','None','faceAlpha',0.2);
    %plot_gaussian_ellipsoid([nanmean(temp(1,dispFlag)) nanmean(temp(2,dispFlag))], cov(temp','omitrows'),1, {'LineWidth',2,'Color',[0.8 0.8 0.8]});

end

tempB = reinfBias(:); tempA = reinfBias(:)-probeBias(:); lm = fitlm(tempB,tempA);
tempX = [0 0.85]; tempY = tempX * lm.Coefficients{2,1} + lm.Coefficients{1,1};

if ~isempty(p.Results.randomComp)
    scatter(tempB, tempA,30,fn_wheelColorsPT('probe'),'filled');
    [r2,p2] = corrcoef(tempB,tempA,'rows','complete');
    %plot([0 0],[-0.8 0.8],'LineWidth',0.5,'Color',[0.8 0.8 0.8]);
    plot([0 1],[0 0],'LineWidth',0.5,'Color',[0.8 0.8 0.8]);
    plot(tempX,tempY,'Color',fn_wheelColorsPT('probe'),'LineWidth',0.5);
    tempAxis = 0:0.01:0.85; [Ypred,YCI] = predict(lm, tempAxis');
    fill([tempAxis fliplr(tempAxis)], [YCI(:,1);flipud(YCI(:,2))]',fn_wheelColorsPT('probe'),'LineStyle','None','faceAlpha',0.2);

    title(['corr = ' num2str(r2(1,2)) ', p = ' num2str(p2(1,2)) ', slope = ' num2str(lm.Coefficients.Estimate(2)) newline...
        'corr = ' num2str(r1(1,2)) ', p = ' num2str(p1(1,2)) ', slope = ' num2str(lm2.Coefficients.Estimate(2))])
    
    %plot_gaussian_ellipsoid([nanmean(tempB) nanmean(tempA)], cov(tempB,tempA,'omitrows'),1, {'LineWidth',2,'Color',matlabColors(1,0.4)});
else
    scatter(tempB, tempA,30,tempColor,'filled');
    plot(tempX,tempY,'Color',[0.6 0.6 0.6],'LineWidth',0.5);
    title(['corr = ' num2str(corr(tempB,tempA,'rows','complete')) ', p = ' num2str(lm.Coefficients.pValue(2))])

end
xlim([0 1]); xticks(0:0.25:1); ylim([-0.8 0.8]);yticks(-0.8:0.4:0.8)




end