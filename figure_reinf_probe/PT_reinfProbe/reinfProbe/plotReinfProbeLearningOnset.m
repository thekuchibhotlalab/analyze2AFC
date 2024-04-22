%function plotReinfProbeLearningOnset(mouseMega)

%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);

%% FIG 1 -- load and align probe data
% align data by minimizing performance variance
probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
outCell = mouseMega.objFun('binProbeByTrialFromLearningOnset',{[reinfThre probeThre],probeTrialBin});
% probe data
probeData = fn_catStructField(2, outCell{1}{:});
probeIdx = cellfun(@fn_getBinMidPoint,outCell{2},'UniformOutput',false);
probeAlignPoint = cell2mat(outCell{4}); 
probeData.probeData = reshape(probeData.probeData,[size(probeData.probeData,1) nMouse 6]);
probeData.befData = reshape(probeData.befData,[size(probeData.probeData,1) nMouse 6]);
probeData.aftData = reshape(probeData.aftData,[size(probeData.probeData,1) nMouse 6]);

[probeDataRawAligned] = getProbeDataCountRaw(outCell{1},probeAlignPoint,-2:1);
[probeDataAlign,probeAlignPoint] = attachNan(probeData, probeAlignPoint);
% reinf data
reinfData = mouseMega.loadReinf;
reinfAlignPoint = cell2mat(outCell{3});
[reinfDataAlign,reinfAlignPoint] = attachNan(reinfData, reinfAlignPoint);
% plot axis
reinfAxis = (1:size(reinfDataAlign.acc,1)) - reinfAlignPoint;
tempProbeBin = probeTrialBin;
probeAxis = (tempProbeBin:tempProbeBin:tempProbeBin*size(probeDataAlign.probeAcc,1)) - tempProbeBin*probeAlignPoint + tempProbeBin/2;

xlimit = [-700 700]; probeSel = probeAlignPoint + (-2:1);

figure; subplot(2,1,1)
plotReinfProbe(reinfDataAlign.acc,probeDataAlign.probeAcc,reinfAxis,probeAxis,probeSel)
plot(xlimit,[0.5 0.5],'LineWidth',2,'Color',[0.8 0.8 0.8])
xlim(xlimit); ylim([0.5 1.0]); yticks(0.5:0.25:1.0); xticks(-600:300:600)

subplot(2,1,2)
plotReinfProbe(reinfDataAlign.bias,abs(probeDataAlign.probeBias),reinfAxis,probeAxis,probeSel)
xlim(xlimit); ylim([0 0.5]); yticks(0:0.1:0.6);xticks(-600:300:600)
%% NEW PLOT 1.1 -- BAR PLOT NEW
probeCount = mouseMega.objFun('computeProbeAlignByTrialBefAft',{'befAft',[reinfThre probeThre],probeTrialBin});
probeCountPerf = mouseMega.objFun('computeProbeAlignByPerfBefAft',{'befAft',[0.6 0.8],400});

%%
dayBef = 2; dayAft = 2;
reinfSelBin = (reinfAlignPoint-dayBef*probeTrialBin):probeTrialBin: (reinfAlignPoint+dayAft*probeTrialBin);
probeSelBin = probeAlignPoint-dayBef: probeAlignPoint+dayAft;

tempAnovaAcc = []; tempAnovaBias = []; anovaLabel1 = []; anovaLabel2 = {};
for i = 1:length(probeSelBin)-1
    %tempAnovaAcc(i,1,:) = nanmean(reinfDataAlign.acc(reinfSelBin(i) : (reinfSelBin(i+1)-1),:),1);
    tempAnovaAcc(i,1,:) = reinfDataAlign.acc( round((reinfSelBin(i)+(reinfSelBin(i+1)-1))/2),:);
    tempAnovaAcc(i,2,:) = probeDataAlign.probeAcc(probeSelBin(i),:);
    %tempAnovaBias(i,1,:) = nanmean(abs(reinfDataAlign.bias(reinfSelBin(i) : (reinfSelBin(i+1)-1),:)),1);
    tempAnovaBias(i,1,:) = abs(reinfDataAlign.bias( round((reinfSelBin(i)+(reinfSelBin(i+1)-1))/2),:));
    tempAnovaBias(i,2,:) = abs(probeDataAlign.probeBias(probeSelBin(i),:));


    anovaLabel1(i,1,:) = ones(1,size(tempAnovaAcc,3))*i;
    anovaLabel1(i,2,:) = ones(1,size(tempAnovaAcc,3))*i;
    tempCell = cell(1,size(tempAnovaAcc,3)); tempCell(:)= {'reinf'};
    anovaLabel2(i,1,:) = tempCell;
    tempCell = cell(1,size(tempAnovaAcc,3)); tempCell(:)= {'probe'};
    anovaLabel2(i,2,:) = tempCell;
end
tempAnovaAcc = tempAnovaAcc(:);
tempAnovaBias = tempAnovaBias(:);
anovaLabel1 = anovaLabel1(:);
anovaLabel2 = anovaLabel2(:);
subject = repmat(reshape(1:13,[1 1 13]),[2 2 1]); subject =subject(:);

%nanFlag = isnan(tempAnovaAcc);
%tempAnovaAcc = tempAnovaAcc(~nanFlag);
%tempAnovaBias = tempAnovaBias(~nanFlag);
%anovaLabel1 = anovaLabel1(~nanFlag);
%anovaLabel2 = anovaLabel2(~nanFlag);
%subject = subject(~nanFlag);

save('C:\Users\zzhu34\Documents\gitRep\analyze2AFCObj\pythonData\reinfProbeData.mat',...
    'anovaLabel1','anovaLabel2','tempAnovaAcc','tempAnovaBias','subject');

[pAcc,tblAcc,statAcc] = anovan(tempAnovaAcc,{anovaLabel1,anovaLabel2});
[pBias,tblBias,statBias] = anovan(tempAnovaBias,{anovaLabel1,anovaLabel2});
[results,~,~,gnames] = multcompare(statBias,'Dimension',[1 2]);
tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A") = gnames(tbl.("Group A"));
tbl.("Group B") = gnames(tbl.("Group B"));
%% FIG 1.1 -- BAR PLOT COMPARISON OF BEFORE VS> AFTER

figure; subplot(2,1,1)
plotBar(probeDataRawAligned.befAcc,probeDataRawAligned.probeAcc,probeDataRawAligned.aftAcc,'left')
ylim([0.5 1.0]); yticks([0.5 0.75 1])
subplot(2,1,2)
plotBar(abs(probeDataRawAligned.befBias),abs(probeDataRawAligned.probeBias),abs(probeDataRawAligned.aftBias),'right')
ylim([0 0.5]); yticks([0 0.25 0.5])

dotColor = [0.7 0.7 0.7];
figure; subplot(2,1,1); hold on;
boxplot([probeDataRawAligned.befAcc,probeDataRawAligned.probeAcc,probeDataRawAligned.aftAcc])
scatter(ones(size(probeDataRawAligned.befAcc))+rand(size(probeDataRawAligned.befAcc))*0.1-0.05,...
        probeDataRawAligned.befAcc,15,dotColor,'filled');
scatter(ones(size(probeDataRawAligned.probeAcc))*2+rand(size(probeDataRawAligned.probeAcc))*0.1-0.05,...
        probeDataRawAligned.probeAcc,15,dotColor,'filled');
scatter(ones(size(probeDataRawAligned.aftAcc))*3+rand(size(probeDataRawAligned.aftAcc))*0.1-0.05,...
        probeDataRawAligned.aftAcc,15,dotColor,'filled');
ylim([0.5 1.0]); yticks([0.5 0.75 1]); xlim([0.5 3.5])
subplot(2,1,2); hold on; 
boxplot(abs([probeDataRawAligned.befBias,probeDataRawAligned.probeBias,probeDataRawAligned.aftBias]))
scatter(ones(size(probeDataRawAligned.befBias))+rand(size(probeDataRawAligned.befBias))*0.2-0.1,...
        abs(probeDataRawAligned.befBias),15,dotColor,'filled');
scatter(ones(size(probeDataRawAligned.probeBias))*2+rand(size(probeDataRawAligned.probeBias))*0.2-0.1,...
        abs(probeDataRawAligned.probeBias),15,dotColor,'filled');
scatter(ones(size(probeDataRawAligned.aftBias))*3+rand(size(probeDataRawAligned.aftBias))*0.2-0.1,...
        abs(probeDataRawAligned.aftBias),15,dotColor,'filled');
ylim([0 0.5]); yticks([0 0.25 0.5]);xlim([0.5 3.5])
%% bar plot by summing across all biased blocks
probeBias = cellfun(@(x)(x.probeBias),probeCount{1}); 
befBias = cellfun(@(x)(x.befBias),probeCount{1}); 
aftBias = cellfun(@(x)(x.aftBias),probeCount{1});  
[a,b,c] = anova1([abs(befBias') abs(probeBias'),abs(aftBias')]);
[a,b,c] = kruskalwallis([abs(befBias') abs(probeBias'),abs(aftBias')]);
d = multcompare(c);

figure;

plotBar(abs(befBias),abs(probeBias),abs(aftBias))
ylim([0 0.5]); yticks([0 0.25 0.5])

%% BAR PLOT COMPARISON BEF VS. AFT
probeByDay = fn_cell2mat(cellfun(@(x)((abs(x.probe))),probeCountPerf{2},'UniformOutput',false),1); 
befByDay = fn_cell2mat(cellfun(@(x)((abs(x.bef))),probeCountPerf{2},'UniformOutput',false),1); 
aftByDay = fn_cell2mat(cellfun(@(x)((abs(x.aft))),probeCountPerf{2},'UniformOutput',false),1); 

probe = fn_cell2mat(cellfun(@(x)(nanmean(abs(x.probe),1)),probeCountPerf{4},'UniformOutput',false),1); 
bef = fn_cell2mat(cellfun(@(x)(nanmean(abs(x.bef),1)),probeCountPerf{4},'UniformOutput',false),1); 
aft = fn_cell2mat(cellfun(@(x)(nanmean(abs(x.aft),1)),probeCountPerf{4},'UniformOutput',false),1); 

figure;
subplot(2,1,2);fn_boxplot({befByDay(:,2)/2+aftByDay(:,2)/2,probeByDay(:,2)},'paired',true,'dotLoc','side',...
    'boxplotArgIn',{}); 
xticks([1 2]); xticklabels({'reinf','probe'}); ylim([0 1.0]); yticks([0 0.25 0.5 0.75 1.0])
[h,p] = signrank(befByDay(:,2)/2+aftByDay(:,2)/2,probeByDay(:,2));  title(['p= ' num2str(h,'%.4f')])
subplot(2,1,1);fn_boxplot({befByDay(:,1)/2+aftByDay(:,1)/2,probeByDay(:,1)},'paired',true,'dotLoc','side',...
    'boxplotArgIn',{}); 
xticks([1 2]); xticklabels({'reinf','probe'}); ylim([0.25 1.0]); yticks([0.25 0.5 0.75 1.0])
[h,p] = signrank(befByDay(:,1)/2+aftByDay(:,1)/2,probeByDay(:,1));  title(['p= ' num2str(h,'%.4f')])

figure;
subplot(2,1,2);fn_boxplot({bef(:,2)/2+aft(:,2)/2,probe(:,2)},'paired',true,'dotLoc','side',...
    'boxplotArgIn',{}); 
xticks([1 2]); xticklabels({'reinf','probe'}); ylim([0 0.6]); yticks([0 0.3 0.6])
[h,p] = signrank(bef(:,2)/2+aft(:,2)/2,probe(:,2));  title(['p= ' num2str(h,'%.4f')])
subplot(2,1,1);fn_boxplot({bef(:,1)/2+aft(:,1)/2,probe(:,1)},'paired',true,'dotLoc','side',...
    'boxplotArgIn',{}); 
xticks([1 2]); xticklabels({'reinf','probe'}); ylim([0.5 1.0]); yticks([0.5 0.75 1.0])
[h,p] = signrank(bef(:,1)/2+aft(:,1)/2,probe(:,1)); title(['p= ' num2str(h,'%.4f')])

figure;fn_boxplot({bef(:,2),probe(:,2),aft(:,2)},'paired',true); xticklabels({'bef','probe','aft'})
[p,tbl,stats] = kruskalwallis([(bef(:,2)) (probe(:,2)) (aft(:,2))]);
d = multcompare(stats);

figure;fn_boxplot({befByDay(:,2),probeByDay(:,2),aftByDay(:,2)},'paired',true); xticklabels({'bef','probe','aft'})
[p,tbl,stats] = kruskalwallis([(befByDay(:,2)) (probeByDay(:,2)) (aftByDay(:,2))]);
d = multcompare(stats);
%% bar plot by subsampling similar amount of trials with probe
probe = fn_cell2mat(cellfun(@(x)(nanmean(abs(x.probe),1)),probeCount{2},'UniformOutput',false),1); 
bef = fn_cell2mat(cellfun(@(x)(nanmean(abs(x.bef),1)),probeCount{2},'UniformOutput',false),1); 
aft = fn_cell2mat(cellfun(@(x)(nanmean(abs(x.aft),1)),probeCount{2},'UniformOutput',false),1); 
[a,b,c] = anova1([(bef(:,2)) (probe(:,2)) (aft(:,2))]);
d = multcompare(c);

figure;
plotBar(bef(:,2),probe(:,2),aft(:,2))
ylim([0 0.5]); yticks([0 0.25 0.5])

probe = fn_cell2mat(cellfun(@(x)(nanmean(abs(x.probe),1)),probeCount{4},'UniformOutput',false),1); 
bef = fn_cell2mat(cellfun(@(x)(nanmean(abs(x.bef),1)),probeCount{4},'UniformOutput',false),1); 
aft = fn_cell2mat(cellfun(@(x)(nanmean(abs(x.aft),1)),probeCount{4},'UniformOutput',false),1); 
[a,b,c] = anova1([(bef(:,2)) (probe(:,2)) (aft(:,2))]);
d = multcompare(c);

[p,tbl,stats] = kruskalwallis([(bef(:,2)) (probe(:,2)) (aft(:,2))]);
c = multcompare(stats);
figure;
plotBar(bef(:,2),probe(:,2),aft(:,2))
ylim([0 0.7]); yticks([0 0.25 0.5])

figure;fn_boxplot({bef(:,2),probe(:,2),aft(:,2)},'paired',true); xticklabels({'bef','probe','aft'})

probeByDay = fn_cell2mat(cellfun(@(x)((abs(x.probe))),probeCount{4},'UniformOutput',false),1); 
befByDay = fn_cell2mat(cellfun(@(x)((abs(x.bef))),probeCount{4},'UniformOutput',false),1); 
aftByDay = fn_cell2mat(cellfun(@(x)((abs(x.aft))),probeCount{4},'UniformOutput',false),1); 

[p,tbl,stats] = kruskalwallis([(befByDay(:,2)) (probeByDay(:,2)) (aftByDay(:,2))]);
c = multcompare(stats)

[p,tbl,stats] = anova1([(befByDay(:,2)) (probeByDay(:,2)) (aftByDay(:,2))]);
c = multcompare(stats)

[h,p] = signrank((probeByDay(:,2)),(befByDay(:,2)+aftByDay(:,2))/2);
[h,p] = kstest(probeByDay(:,2));
figure; fn_plotBarPaired({befByDay(:,2)/2+aftByDay(:,2)/2,probeByDay(:,2)})
figure; plotBar(befByDay(:,2),probeByDay(:,2),aftByDay(:,2))

figure;fn_boxplot({befByDay(:,2),probeByDay(:,2),aftByDay(:,2)},'paired',true); xticklabels({'bef','probe','aft'})




figure; hold on;
plot([0 0],[-0.6 0.6],'Color',[0.8 0.8 0.8]); plot([-0.8 0.8],[0 0],'Color',[0.8 0.8 0.8]);
scatter((befByDay(:,2)+aftByDay(:,2))/2-(probeByDay(:,2)),-(befByDay(:,1)+aftByDay(:,1))/2+(probeByDay(:,1)),...
    20,matlabColors(1),'filled'); 
xlim([-0.8 0.8]);ylim([-0.6 0.6])

fitlm((befByDay(:,2)+aftByDay(:,2))/2-(probeByDay(:,2)),-(befByDay(:,1)+aftByDay(:,1))/2+(probeByDay(:,1)))

fitlm((befByDay(:,2)+aftByDay(:,2))/2,(befByDay(:,2)+aftByDay(:,2))/2-(probeByDay(:,2)))
%%
[a,b,c] = anova1([abs(probeDataRawAligned.befAcc) abs(probeDataRawAligned.probeAcc),abs(probeDataRawAligned.aftAcc)]);
d = multcompare(c);
%[h,p] = kstest(abs(probeDataRawAligned.befBias))
%[h,p] = kstest(abs(probeDataRawAligned.probeBias))
%[h,p] = kstest(abs(probeDataRawAligned.aftBias))

[a,b,c] = anova1([abs(probeDataRawAligned.befBias) abs(probeDataRawAligned.probeBias),abs(probeDataRawAligned.aftBias)]);
d = multcompare(c);

temp = [nanmean(abs(probeDataAlign.befBias(probeSel,:)),1)' nanmean(abs(probeDataAlign.probeBias(probeSel,:)),1)' ...
    nanmean(abs(probeDataAlign.aftBias(probeSel,:)),1)'];
[a,b,c] = anova1([nanmean(abs(probeDataAlign.befBias(probeSel,:)),1)' nanmean(abs(probeDataAlign.probeBias(probeSel,:)),1)' ...
    nanmean(abs(probeDataAlign.aftBias(probeSel,:)),1)']);
d = multcompare(c);
figure;
fn_plotBarPaired({abs(probeDataRawAligned.befAftBias),abs(probeDataRawAligned.probeBias)})
ylim([0 0.5]); yticks([0 0.25 0.5])
%% FIG 1 SUPPLEMENT -- FIRST VS. LAST HALF OF PROBE DATA
figure; subplot(2,2,1)
plotReinfProbe(reinfDataAlign.acc,probeDataAlign.firstHalfAcc,reinfAxis,probeAxis)
xlim(xlimit); ylim([0.4 0.9]); title('First Half Probe')

subplot(2,2,2)
plotReinfProbe(reinfDataAlign.acc,probeDataAlign.lastHalfAcc,reinfAxis,probeAxis)
xlim(xlimit); ylim([0.4 0.9]); title('Last Half Probe')

subplot(2,2,3)
plotReinfProbe(reinfDataAlign.bias,abs(probeDataAlign.firstHalfBias),reinfAxis,probeAxis)
xlim(xlimit); ylim([0 0.6]); 

subplot(2,2,4)
plotReinfProbe(reinfDataAlign.bias,abs(probeDataAlign.lastHalfBias),reinfAxis,probeAxis)
xlim(xlimit); ylim([0 0.6])

%% FIG 1 SUPPLEMENT -- dPRIME
figure;
plotReinfProbe(reinfDataAlign.dp,probeDataAlign.probeDp,reinfAxis,probeAxis)
xlim(xlimit); ylim([-0.5 3])
%% FIG 1 SUPPLEMENT -- ALIGNMENT
figure; subplot(2,1,2); hold on;
for i = 1:nMouse; plot(reinfAxis,reinfDataAlign.acc(:,i),'Color',fn_wheelColorsPT('Reinf',0.5),'LineWidth',0.8); end 
plot(reinfAxis,nanmean(reinfDataAlign.acc,2),'Color',fn_wheelColorsPT('Reinf'),'LineWidth',2);
xlim([-1000 2000])

subplot(2,1,1); hold on;
tempAcc = mouseMega.getProp('behav','field','acc','matFlag',true);
for i = 1:nMouse; plot(tempAcc(:,i),'Color',[0.6 0.6 0.6],'LineWidth',0.8); end 
plot(nanmean(tempAcc,2),'Color',[0 0 0],'LineWidth',1.5);
xlim([0 3000])

figure;
subplot(2,3,5);
imagesc(reinfDataAlign.acc');tempZero = find(reinfAxis==0); xticks([tempZero tempZero+2000]); xticklabels({'0','2000'}); xlim([-1000 2500]+tempZero)

subplot(2,3,3); hold on;
for i = 1:nMouse; plot(probeAxis,probeDataAlign.probeAcc(:,i),'Color',fn_wheelColorsPT('Probe',0.5),'LineWidth',0.8); end 
plot(probeAxis,nanmean(probeDataAlign.probeAcc,2),'Color',fn_wheelColorsPT('Probe'),'LineWidth',2);
xlim([-1000 2500])

subplot(2,3,6); hold on;
for i = 1:nMouse; plot(reinfAxis,reinfDataAlign.acc(:,i),'Color',fn_wheelColorsPT('Reinf',0.4),'LineWidth',0.8); end 
for i = 1:nMouse; plot(probeAxis,probeDataAlign.probeAcc(:,i),'Color',fn_wheelColorsPT('Probe',0.8),'LineWidth',0.8); end 
xlim([-1000 2500])

subplot(2,3,2);
imagesc(tempAcc'); 


selAni1 = 8; selAni2 = 11;
figure; subplot(2,1,1); 
plot([(1608-150) (5557+1608-150-1)],[0.7 0.7],'Color',[0.8 0.8 0.8],'LineWidth',2); hold on;
idx1 = find(reinfData.acc(:,selAni1)>=0.7,1); idx2 = find(reinfData.acc(:,selAni2)>=0.7,1);
plot([(1608-150)+idx1 (1608-150)+idx1],[0.4 1],'-','Color',[matlabColors(1,0.4) 0.8],'LineWidth',2)
plot([(1608-150)+idx2 (1608-150)+idx2],[0.4 1],'-','Color',[matlabColors(2,0.4) 0.8],'LineWidth',2)

plot((1608-150):(5557+1608-150-1),reinfData.acc(:,selAni1),'LineWidth',2, 'Color',matlabColors(1)); 
plot((1608-150):(5557+1608-150-1),reinfData.acc(:,selAni2),'LineWidth',2,'Color',matlabColors(2))

xlim([1500 4200]); xticks([1500 4000]); xticklabels({'0','2500'})
ylim([0.4 1]); yticks(0.4:0.3:1)

subplot(2,1,2);
plot([(1608-150) (5557+1608-150-1)],[0.7 0.7],'Color',[0.8 0.8 0.8],'LineWidth',2); hold on;
idx1 = find(reinfDataAlign.acc(:,selAni1)>=0.7,1); idx2 = find(reinfDataAlign.acc(:,selAni2)>=0.7,1);
plot([idx1 idx1],[0.4 1],'-','Color',[matlabColors(1,0.4) 0.8],'LineWidth',2)
plot([idx2 idx2],[0.4 1],'-','Color',[matlabColors(2,0.4) 0.8],'LineWidth',2)

plot(reinfDataAlign.acc(:,selAni1),'LineWidth',2,'Color',matlabColors(1)); 
plot(reinfDataAlign.acc(:,selAni2),'LineWidth',2,'Color',matlabColors(2))
xlim([1500 4200]); xticks([idx1-1000 idx1 idx1+1000]); xticklabels({'-1000','0','1000'})
ylim([0.4 1]); yticks(0.4:0.3:1)

figure; %subplot(1,2,1); plot(sum(isnan(reinfDataAlign.acc),2))
a = sum(isnan(reinfDataAlign.acc),2); varianceStart = find(a==5,1,'first'); varianceEnd = find(a==5,1,'last');
%varianceStart = 1828; varianceEnd = 4926;
varAlign = nanstd(reinfDataAlign.acc(varianceStart:varianceEnd,:),0,2);
varNoAlign = nanstd(reinfData.acc(1:3000,:),0,2);

varAlign = [varAlign;nan(16,1)];
varAlign = fn_binMean(varAlign,50);
varNoAlign = fn_binMean(varNoAlign,50);
%subplot(1,2,2); hold on;
p = fn_plotComparison({varAlign,varNoAlign},'compType','errorbarWithDot','paired',false,'dotType','random');
%bar([nanmean(varAlign); nanmean(varNoAlign)],'EdgeColor',[0 0 0],'FaceColor','None'); 
%errorbar([1 2], [nanmean(varAlign), nanmean(varNoAlign)],...
%    [nanstd(varAlign)/sqrt(length(varAlign)), nanstd(varNoAlign)/length(varNoAlign)],...
%    'LineStyle','None','Color',[0 0 0]);
ylim([0 0.15]);xlim([0.6 2.4]); xticks([1 2]); xticklabels({'aligned','not aligned'})
%% FIG 1 SUPPLEMENT -- INDIVIDUAL ANIMALS
for i = 1:nMouse
    maxidx = find(~isnan(reinfDataAlign.acc(:,i)),1,'last');
    minidx = find(~isnan(reinfDataAlign.acc(:,i)),1,'first');
    figure; subplot(2,1,1);
    plotReinfProbeIndiv(reinfDataAlign.acc(:,i),probeDataAlign.probeAcc(:,i),reinfAxis,probeAxis)
    %plot([reinfAxis(1) reinfAxis(end)],[0.5 0.5],'LineWidth',2,'Color',[0.8 0.8 0.8])
    xlim([reinfAxis(minidx) reinfAxis(maxidx)]); ylim([0.4 1.0])
    %subplot(3,1,2);
    %plotReinfProbeIndiv(reinfDataAlign.dp(:,i),probeDataAlign.probeDp(:,i),reinfAxis,probeAxis)
    %xlim([reinfAxis(minidx) reinfAxis(maxidx)]); ylim([0 3])
    subplot(2,1,2)
    plotReinfProbeIndiv(reinfDataAlign.biasDir(:,i),(probeDataAlign.probeBias(:,i)),reinfAxis,probeAxis)
    xlim([reinfAxis(minidx) reinfAxis(maxidx)]); ylim([-1 1.0])
    
end

%% FIG 1 SUPPLEMENT -- 2 INDIVIDUAL ANIMALS EXAMPLES
for mouse1 = 1:13
    maxidx = find(~isnan(reinfDataAlign.acc(:,mouse1)),1,'last');
    minidx = find(~isnan(reinfDataAlign.acc(:,mouse1)),1,'first');
    
    minidx = max([minidx reinfAlignPoint-700]);
    maxidx = min([maxidx reinfAlignPoint+700]);
    
     probeSel = probeAlignPoint + (-2:1);
    
    figure; subplot(2,1,1);
    plotReinfProbeIndiv(smoothdata(reinfDataAlign.acc(:,mouse1),'gaussian',30),...
        probeDataAlign.probeAcc(:,mouse1),reinfAxis,probeAxis,probeSel)
    xlim([reinfAxis(minidx) reinfAxis(maxidx)]); ylim([0.2 1.0])
    
    subplot(2,1,2)
    hold on; plot([reinfAxis(minidx), reinfAxis(maxidx)],[0 0],'Color',[0.8 0.8 0.8],'LineWidth',2)
    plotReinfProbeIndiv(smoothdata(abs(reinfDataAlign.biasDir(:,mouse1)),'gaussian',30),abs(probeDataAlign.probeBias(:,mouse1)),reinfAxis,probeAxis,probeSel)
    xlim([reinfAxis(minidx) reinfAxis(maxidx)]); ylim([0 1.0])
    exportgraphics(gcf,'D:\OneDrive - Johns Hopkins University\Lab meeting\reinf_probe_paper_fig\reinf_probe\probe_fig.pdf','Append',true);
end

%% FIG 2 -- CORRELATION BETWEEN BIAS AND PERFORMANCE
probeSelIdx = probeAlignPoint + (-2:1);
probeAccSel = (probeDataAlign.probeAcc(probeSelIdx,:)); 
befAccSel = (probeDataAlign.befAcc(probeSelIdx,:));
aftAccSel = (probeDataAlign.aftAcc(probeSelIdx,:)); 

probeBiasSel = (abs(probeDataAlign.probeBias(probeSelIdx,:))); 
befBiasSel = (abs(probeDataAlign.befBias(probeSelIdx,:)));
aftBiasSel = (abs(probeDataAlign.aftBias(probeSelIdx,:))); 

tempAcc = []; tempAcc(1,:,:) = (befAccSel+aftAccSel)/2; tempAcc(2,:,:) = probeAccSel;
tempBias = []; tempBias(1,:,:) = (befBiasSel+aftBiasSel)/2; tempBias(2,:,:) = probeBiasSel;
day = repmat(1:4,[2 1 13]); context = []; context(1,:,:) = zeros(4,13); context(2,:,:) = ones(4,13);
subject = repmat(reshape(1:13,[1 1 13]),[2 4 1]);
tempAcc = tempAcc(:); tempBias = tempBias(:); day = day(:); subject = subject(:);  context = context(:);
save('C:\Users\zzhu34\Documents\gitRep\analyze2AFCObj\pythonData\reinfProbeData_befAft.mat', ...
    'tempAcc','tempBias','day','subject','context');

colors = parula; tempIdx = round(linspace(1, size(colors,1),nMouse+1));
colors = colors(tempIdx(1:end-1),:);
tempColor = []; tempColor(1,:,:) = colors; tempColor = repmat(tempColor,size(probeAccSel,1),1,1);
tempColor = reshape(tempColor,size(tempColor,1)* size(tempColor,2),[]);

a = probeAccSel(:)- befAccSel(:)/2 - aftAccSel(:)/2;
b = befBiasSel(:)/2 + aftBiasSel(:)/2 - probeBiasSel(:);
d = befBiasSel(:)/2 + aftBiasSel(:)/2;

figure; 
subplot(3,1,1); hold on;
scatter(d,a,50,tempColor,'filled'); c = fitlm(d,a);  xlabel('bias Reinf'); ylabel('acc Probe - Reinf')
title(['accP-R vs. biasR, p= ' num2str(c.Coefficients.pValue(2),'%.4f')])
xlim([0 0.8]); ylim([-0.5 0.5])

subplot(3,1,2); hold on;
c = fitlm(d,b); 
tempX = 0:0.01:1.0; tempY = tempX * c.Coefficients{2,1} + c.Coefficients{1,1};
plot(tempX,tempY,'Color',[0.8 0.8 0.8])
scatter(d,b,50,tempColor,'filled'); xlabel('bias Reinf'); ylabel('bias Probe - Reinf')
title(['biasR-P vs. biasR, p= ' num2str(c.Coefficients.pValue(2),'%.4f')])
xlim([0 1.0]); ylim([-0.8 0.8]); yticks(-0.8:0.4:0.8);xticks(0:0.2:1);

subplot(3,1,3); hold on;
c = fitlm(b,a);
tempX = -0.8:0.01:0.8; tempY = tempX * c.Coefficients{2,1} + c.Coefficients{1,1};
plot(tempX,tempY,'Color',[0.8 0.8 0.8])
scatter(b,a,50,tempColor,'filled');  xlabel('bias Probe-Reinf'); ylabel('acc Probe - Reinf')
title(['accP-R vs. biasR-P, p= ' num2str(c.Coefficients.pValue(2),'%.4f')])
xlim([-0.8 0.8]); ylim([-0.4 0.4]); xticks(-0.8:0.4:0.8); yticks(-0.8:0.4:0.8)

figure; [h,p] = ttest(befAccSel(:)/2 + aftAccSel(:)/2,probeAccSel(:),'tail','left');
subplot(2,1,1); hold on;plot([0.4 1.0],[0.4 1.0],'Color',[0.8 0.8 0.8])
scatter(befAccSel(:)/2 + aftAccSel(:)/2,probeAccSel(:),50,tempColor,'filled'); xlabel('acc Reinf'); ylabel('acc Probe')
title(['accR vs. accP, p = ' num2str(p,'%.4f')])
xlim([0.4 1.0]); ylim([0.4 1.0]); 

[h,p] = ttest(befBiasSel(:)/2 + aftBiasSel(:)/2,probeBiasSel(:),'tail','right');
subplot(2,1,2); hold on;plot([0 0.8],[0 0.8],'Color',[0.8 0.8 0.8])
scatter(befBiasSel(:)/2 + aftBiasSel(:)/2,probeBiasSel(:),50,tempColor,'filled'); xlabel('bias Reinf'); ylabel('bias Probe')
title(['biasR vs. biasP, p =' num2str(p,'%.4f')])
xlim([0 0.8]); ylim([0 0.8])
%% FIG 2 SUPPLEMENT -- INDIVIDUAL ANIMAL EXAMPLE
selMouse = 5; tempIdx = 4*(selMouse-1)+1:4*selMouse;
figure;
subplot(3,1,1); hold on;plot([0.4 1.0],[0.4 1.0],'Color',[0.8 0.8 0.8])
scatter(befAccSel(:,selMouse)/2 + aftAccSel(:,selMouse)/2,probeAccSel(:,selMouse),50,tempColor(tempIdx,:),'filled'); xlabel('acc Reinf'); ylabel('acc Probe')
title(['accR vs. accP, mouse ' int2str(selMouse)])
xlim([0.4 1.0]); ylim([0.4 1.0])

subplot(3,1,2); hold on;plot([0 0.8],[0 0.8],'Color',[0.8 0.8 0.8])
scatter(befBiasSel(:,selMouse)/2 + aftBiasSel(:,selMouse)/2,probeBiasSel(:,selMouse),50,tempColor(tempIdx,:),'filled'); xlabel('bias Reinf'); ylabel('bias Probe')
title(['biasR vs. biasP, mouse ' int2str(selMouse)])
xlim([0 0.8]); ylim([0 0.8])

subplot(3,1,3); scatter(ones(1,nMouse),nMouse:-1:1,50,colors,'filled');
%% FIG 2 SUPPLEMENT -- ERROBAR OF AVERAGED by animal
probeSelIdx = probeAlignPoint+(-1:0);

probeAccSel = probeDataAlign.probeAcc(probeSelIdx,:); 
befAccSel = probeDataAlign.befAcc(probeSelIdx,:);
aftAccSel = probeDataAlign.aftAcc(probeSelIdx,:); 

probeBiasSel = abs(probeDataAlign.probeBias(probeSelIdx,:)); 
befBiasSel = abs(probeDataAlign.befBias(probeSelIdx,:));
aftBiasSel = abs(probeDataAlign.aftBias(probeSelIdx,:)); 

nData = sum(~isnan(probeAccSel),1);


figure;
[~,p] = ttest( nanmean(befAccSel,1)/2 + nanmean(aftAccSel,1)/2,nanmean(probeAccSel,1),'tail','left');
subplot(2,1,1); hold on;
plotScatterErrorBar(befAccSel/2 + aftAccSel/2,probeAccSel);
plot([0.4 1.0],[0.4 1.0],'Color',[0.8 0.8 0.8])
%scatter(nanmean(befAccSel,1)/2 + nanmean(aftAccSel,1)/2,nanmean(probeAccSel,1),50,'filled'); xlabel('acc Reinf'); ylabel('acc Probe')
title(['accR vs. accP, p = ' num2str(p,'%.4f')])
xlim([0.4 1.0]); ylim([0.4 1.0])

[h,p] = ttest(nanmean(befBiasSel,1)/2 + nanmean(aftBiasSel,1)/2,nanmean(probeBiasSel,1),'tail','right');
subplot(2,1,2); hold on;
plotScatterErrorBar(befBiasSel/2 + aftBiasSel/2,probeBiasSel);
plot([0 0.8],[0 0.8],'Color',[0.8 0.8 0.8])
xlabel('bias Reinf'); ylabel('bias Probe')
title(['biasR vs. bias, p =' num2str(p,'%.4f')])
xlim([0 0.8]); ylim([0 0.8])

a = probeAccSel- befAccSel/2 - aftAccSel/2;
b = probeBiasSel- befBiasSel/2 - aftBiasSel/2;
d = befBiasSel/2 + aftBiasSel/2;

figure;
subplot(3,1,1)
p = plotScatterErrorBar(d,a);
xlabel('bias Reinf'); ylabel('acc Probe - Reinf')
title(['accP-R vs. biasR, p= ' num2str(p,'%.4f')])
xlim([0 0.7]); ylim([-0.3 0.3])

subplot(3,1,2)
p = plotScatterErrorBar(d,b);
xlabel('bias Reinf'); ylabel('bias Probe - Reinf')
title(['biasP-R vs. biasR, p= '  num2str(p,'%.4f')])
xlim([0 0.7]); ylim([-0.5 0.5])

subplot(3,1,3)
p = plotScatterErrorBar(b,a);
xlabel('bias Probe-Reinf'); ylabel('acc Probe - Reinf')
title(['accP-R vs. biasP-R, p= '  num2str(p,'%.4f')])
xlim([-0.4 0.4]); ylim([-0.3 0.3])
%% FIG 3 -- ACTION BIAS AND LEARNING RATE CORRELATION
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice); tempOps.learningCurveBin = 100;
allMouse = fn_getObjPT(mice,tempOps);
mouseMega = wheel2AFCmega(allMouse);
% align data by minimizing performance variance
probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
outCell = mouseMega.objFun('binProbeByTrialFromLearningOnset',{[reinfThre probeThre],probeTrialBin});
%%
accThre = 0.80; tempAccIdx = []; tempBias = [];tempLearningOnset = -250;
reinfData = mouseMega.loadReinf; 
reinfData.context = mouseMega.getProp('behav','field','context','matFlag',true);
reinfData.stimProb = mouseMega.getProp('behav','field','stimulus','matFlag',true);
reinfData.stimProb = smoothdata(reinfData.stimProb==1,1,'movmean',100);
tempReinfAlignPoint = cell2mat(outCell{3});
[reinfDataAlign,reinfAlignPoint] = attachNan(reinfData, tempReinfAlignPoint);
for i = 1:size(reinfData.acc,2)
    tempAccIdx(i) = find(reinfDataAlign.acc(:,i)>accThre,1,'first')-reinfAlignPoint;
    tempBias(i) = nanmean(abs(reinfDataAlign.bias(reinfAlignPoint+tempLearningOnset:reinfAlignPoint+tempAccIdx(i),i)));
    %disp(accIdx(i))
end
disp(tempAccIdx)
figure; c = fitlm(tempBias, tempAccIdx-tempLearningOnset);
scatter(tempBias, tempAccIdx-tempLearningOnset)

%%
colors = parula; tempIdx = round(linspace(1, size(colors,1),nMouse+1));
colors = colors(tempIdx(1:end-1),:);

tempBias = []; tempCorr = []; tempStim = [];
tempLearningOnset = -250;
tempThreshold = 0.75; % use 'double smoothing here', which is not good, change later
figure;
plot(reinfDataAlign.acc)

for i = 1:size(reinfDataAlign.acc,2)
    
    try
        %figure; 
        [param,stat]=sigm_fit(1:size(reinfDataAlign.acc,1),reinfDataAlign.acc(:,i),[],[0.5 0.8 2000 0],0);
        xFit = reinfAlignPoint-1000: reinfAlignPoint+4000;
        yFit = param(1)+(param(2)-param(1))./(1+10.^((param(3)-xFit)*param(4)));
        threTop = find(yFit>0.75,1,'first'); threBot = find(yFit>0.55,1,'first');
        trialToLearn(i) = threTop-threBot;
        tempBias(i) = nanmean(abs(reinfDataAlign.bias(reinfAlignPoint-1000+(threBot:threTop),i)));
        tempCorr(i) = nanmean(reinfDataAlign.context(reinfAlignPoint-1000+(threBot:threTop),i)==2);
        tempStim(i) = nanmean(abs(reinfDataAlign.stimProb(reinfAlignPoint-1000+(threBot:threTop),i)-0.5));
    catch
    end
end

figure; c = fitlm(tempBias([1:4 6:12]), trialToLearn([1:4 6:12]));
subplot(1,3,1);scatter(tempBias(1:12), trialToLearn, 20, colors(1:12,:),'filled')
subplot(1,3,2);scatter(tempCorr(1:12), trialToLearn, 20, colors(1:12,:),'filled')
subplot(1,3,3);scatter(tempStim(1:12), trialToLearn, 20, colors(1:12,:),'filled')

%%
exampleAnimal = 8; tempProbeOn = 796; probeBefAft = 30;
plotTrial = tempProbeOn-probeBefAft : tempProbeOn + probeBefAft + 9; 
figure; subplot(1,2,2)
fn_plotWheelByTrial(mouseMega,exampleAnimal,plotTrial,'probeIdx',probeBefAft+1:probeBefAft+10);

subplot(1,2,1); 
fn_plotChoiceByTrialVertical(mouseMega.mouseCell{exampleAnimal},plotTrial,'csize',15,'probeIdx',probeBefAft+1:probeBefAft+10);
xlim([0.5 2.5])

%%
exampleAnimal = 12; tempProbeOn = 1612; probeBefAft = 15;
plotTrial = tempProbeOn-probeBefAft : tempProbeOn + probeBefAft + 9; 
figure; subplot(2,2,2)
fn_plotWheelByTrial(mouseMega,exampleAnimal,plotTrial,'probeIdx',probeBefAft+1:probeBefAft+10);

subplot(2,2,1); 
fn_plotChoiceByTrialVertical(mouseMega.mouseCell{exampleAnimal},plotTrial,'csize',15,'probeIdx',probeBefAft+1:probeBefAft+10);
xlim([0.5 2.5])

tempAcc = [0.75 1 0.665]; tempBias = [0.5 0 0.67]; 

subplot(2,2,3); hold on; 
plot(1:3,tempAcc,'Color',[0.6 0.6 0.6],'LineWidth',1);  
scatter(1:3,tempAcc,20,[0.6 0.6 0.6],'filled'); xlim([0.5 3.5]); ylim([0.5 1]); yticks([0.5 0.75 1])

subplot(2,2,4); hold on; 
plot(1:3,tempBias,'Color',[0.6 0.6 0.6],'LineWidth',1); 
scatter(1:3,tempBias,20,[0.6 0.6 0.6],'filled'); xlim([0.5 3.5]); ylim([0 1]); yticks([0 0.5 1])
%%
function [probeData,maxAlignPoint] = attachNan(probeData, alignPoint)
    attachDim = 1; 
    maxAlignPoint = max(alignPoint);
    tempFieldNames = fieldnames(probeData);
    for i = 1:length(tempFieldNames)
        tempField = probeData.(tempFieldNames{i});
        if isnumeric(tempField)
            tempMatSize = size(tempField);
            tempMatSize(attachDim) = tempMatSize(attachDim) + maxAlignPoint;
            tempMat = nan(tempMatSize);
            for j = 1:length(alignPoint)
                startPoint = maxAlignPoint-alignPoint(j)+1;
                tempMat(startPoint:startPoint+size(tempField,attachDim)-1,j) = tempField(:,j);

            end
            probeData.(tempFieldNames{i}) = tempMat;
        end        
    end
end

function plotReinfProbe(reinf,probe,reinfAxis,probeAxis,probeSel)
hold on;  
reinfAcc_mean = nanmean(reinf,2); reinfAcc_sem = fn_nansem(reinf,2);
probeAcc_mean = nanmean(probe(probeSel,:),2); probeAcc_sem = fn_nansem(probe(probeSel,:),2);
plot([0 0],[-5 5],'--','Color',[0.8 0.8 0.8],'LineWidth',2);
plot(reinfAxis,reinfAcc_mean,'Color',fn_wheelColorsPT('Reinf'),'LineWidth',2);
f_errorbar = fn_plotFillErrorbar(reinfAxis,reinfAcc_mean,reinfAcc_sem,...
        fn_wheelColorsPT('Reinf'),'faceAlpha',0.3,'LineStyle','none');

plot(probeAxis(probeSel),probeAcc_mean,'--o','Color',fn_wheelColorsPT('Probe'),'LineWidth',2);
f_errorbar = fn_plotFillErrorbar(probeAxis(probeSel),probeAcc_mean,probeAcc_sem,...
        fn_wheelColorsPT('Probe'),'faceAlpha',0.3,'LineStyle','none'); 
end

function plotReinfProbeIndiv(reinf,probe,reinfAxis,probeAxis,probeSel)
hold on;  
plot([0 0],[-5 5],'--','Color',[0.8 0.8 0.8],'LineWidth',2);
plot(reinfAxis,reinf,'Color',fn_wheelColorsPT('Reinf'),'LineWidth',2);
plot(probeAxis(probeSel),probe(probeSel),'--o','Color',fn_wheelColorsPT('Probe'),'LineWidth',2);

end

function p = plotScatterErrorBar(x,y)
meanX = nanmean(x,1); meanY = nanmean(y,1);
errorX = fn_nansem(x,1); errorY = fn_nansem(y,1);
hold on; scatter(meanX,meanY,50,'filled');
errorbar(meanX,meanY,errorX,'LineStyle','none','Color',matlabColors(1));
errorbar(meanX,meanY,errorY,'horizontal','LineStyle','none','Color',matlabColors(1));
c = fitlm(meanX,meanY);
p = c.Coefficients.pValue(2);
    
end

function plotBar(tempBefFlat,tempProbeFlat,tempAftFlat,tail)
    
    nanFlag = isnan(tempProbeFlat); 
    tempBefFlat(nanFlag) = []; tempAftFlat(nanFlag) = []; tempProbeFlat(nanFlag) = [];
    %[hBef,pBef] = ttest(tempBefFlat,tempProbeFlat);
    %[hAft,pAft] = ttest(tempAftFlat,tempProbeFlat);
    [pBef,hBef] = signrank(tempBefFlat,tempProbeFlat);
    [pAft,hAft] = signrank(tempAftFlat,tempProbeFlat);
    bar([nanmean(tempBefFlat) nanmean(tempProbeFlat) nanmean(tempAftFlat)] ,'EdgeColor',[0 0 0],'FaceColor','None'); hold on;
    for i = 1:length(tempBefFlat)
        f = plot([1 2 3],[tempBefFlat(i) tempProbeFlat(i) tempAftFlat(i)],'Color',[0.6 0.6 0.6],'Marker','.','MarkerSize',8,...
            'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',0.5);
    end
    legend(f,['pBef = ' num2str(pBef,'%.2e') newline 'pAft = ' num2str(pAft,'%.2e')],'Location','Best')
    xticks([1 2 3]); xticklabels({'Bef','Probe','Aft'}); xlim([0 4]);
    
end

function [probeDataRawAligned] = getProbeDataCountRaw(probeData,alignPoint,sumIdx)

for i = 1:length(probeData)
    tempIdx = alignPoint(i) + sumIdx; tempIdx(tempIdx<=0) = [];
    tempLen = size(probeData{i}.probeData,1); tempIdx(tempIdx>tempLen) = [];

    tempData.probeData(i,:) = squeeze(nansum(probeData{i}.probeData(tempIdx,:),1));
    tempData.befData(i,:) = squeeze(nansum(probeData{i}.befData(tempIdx,:),1));
    tempData.aftData(i,:) = squeeze(nansum(probeData{i}.aftData(tempIdx,:,:),1));
    tempData.befAftData(i,:) = squeeze(nansum(probeData{i}.befAftData(tempIdx,:),1));

end

[tempData.probeAcc,tempData.probeBias,tempData.probeDP] = fn_getAccBias(tempData.probeData); 
[tempData.befAcc,tempData.befBias,tempData.befDP] = fn_getAccBias(tempData.befData); 
[tempData.aftAcc,tempData.aftBias,tempData.aftDP] = fn_getAccBias(tempData.aftData); 
[tempData.befAftAcc,tempData.befAftBias,tempData.befAftDP] = fn_getAccBias(tempData.befAftData); 
probeDataRawAligned = tempData;


function [acc,bias,dprime] = fn_getAccBias(mat)
    acc1 = mat(:,1)./sum(mat(:,1:2),2); acc2 = mat(:,4)./sum(mat(:,4:5),2);
    acc = (acc1+acc2)/2; bias = acc1-acc2;
    
    dprimeAccThreHigh = 0.95; dprimeAccThreLow = 0.05;
    acc1(acc1>dprimeAccThreHigh) = dprimeAccThreHigh; acc2(acc2>dprimeAccThreHigh) = dprimeAccThreHigh;
    acc1(acc1<dprimeAccThreLow) = dprimeAccThreLow; acc2(acc2<dprimeAccThreLow) = dprimeAccThreLow;
    dprime = norminv(acc1,0,1) - norminv(1-acc2,0,1);
end



end 