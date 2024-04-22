clear;
loadPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\psyTrackFit\';
nPrev = 1;

mouse = {'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
%modelComparison = {'S','SA','SRp','SRn','SRpRn','SARpRn','SB','SBA','SBRp','SBRn','SBRpRn'};
modelComparison = {'S','SA','SRp','SRn','SB'};

for i = 1:length(modelComparison)
    for j = 1:length(mouse)
        if strcmp(modelComparison{i},'S') || strcmp(modelComparison{i},'SB')
            load([loadPath filesep mouse{j} 'psytrack_' modelComparison{i} '.mat'])
        else
            load([loadPath filesep mouse{j} 'psytrack_' modelComparison{i} '_nPrev' int2str(nPrev) '.mat'])
        end
        loglikeli(i,j) = xval_logli/length(xval_pL);
    end
end

loglikeli = loglikeli(2:end,:) - repmat(loglikeli(1,:),size(loglikeli,1)-1,1);
figure;
bar(mean(loglikeli,2),'EdgeColor',[1 1 1],'FaceColor',matlabColors(2,0.9)); hold on;

for i = 1:size(loglikeli,1)
    scatter(i*ones(size(loglikeli,2),1) ,loglikeli(i,:),10,[0 0 0],'filled');
end

xticklabels(modelComparison(2:end))
xtickangle(45); ylabel('log-likeli increase over stimulus-only model')
title(['Prev History ' int2str(nPrev)]); ylim([0 0.12])

%% -------------- PLOT FOR COSYNE-----------------------
clear;
loadPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\psyTrackFit\';
nPrev = 1;
mouse = {'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
%mouse = {'zz054','zz062','zz063','zz066','zz067','zz068','zz069'};
modelComparison = {'S','SARpRn','SARpRn','SB','SBARpRn', 'SBARpRn'};
params = [1 1+3*nPrev 1+3*3 1+1 1+3*nPrev+1 1+3*3+1];

loglikeli = [];
for i = 1:length(modelComparison)
    for j = 1:length(mouse)
        if strcmp(modelComparison{i},'S') || strcmp(modelComparison{i},'SB')
            load([loadPath filesep mouse{j} 'psytrack_' modelComparison{i} '.mat'])
        elseif i==6 || i==3
            load([loadPath filesep mouse{j} 'psytrack_' modelComparison{i} '_nPrev' int2str(3) '.mat'])         
        else
            load([loadPath filesep mouse{j} 'psytrack_' modelComparison{i} '_nPrev' int2str(nPrev) '.mat'])
        end
        loglikeli(i,j) = (-xval_logli)/length(xval_pL);
        BIC(i,j) = -xval_logli * 2 + log(length(xval_logli)) * params(i);
    end
end

BIC = BIC - repmat(BIC(1,:),[size(BIC,1) 1]); BIC = BIC(2:end,:);
fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.3); hold on;
fn_plotComparison(-BIC,'paired',false, 'dotType', 'random', 'compType','errorbarWithDot','barType','none','errorbarArgIn',...
    {'Color',[0.2 0.2 0.2],'LineWidth',1.5,'LineStyle','none'},'scatterArgIn',{15, [0.6 0.6 0.6]},...
    'barplotArgIn',{0.4,'EdgeColor','none','FaceColor',matlabColors(2),'LineStyle','none'});
plot([0.5 5.5],[0 0],'LineWidth',1.5,'Color',[0.8 0.8 0.8]);
xlim([0.5 5.5]); ylim([-100 800]); yticks(0:200:800); xticks(1:5)