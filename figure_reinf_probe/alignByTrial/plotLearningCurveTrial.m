function plotLearningCurveTrial(mouse,selectProtocol)
if nargin == 1; selectProtocol = {'puretone'}; end
global sep;
sep = '\';
rootPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\';
loadPath = [rootPath sep 'trialData\'];
figPath = [rootPath sep 'figure\' sep 'learningCurveTrial' sep]; mkdir(figPath);
trialData = fn_selectProtocol(loadPath, mouse, selectProtocol);

outmat = alignDaysKeepMiss(trialData);

f = plotAnimalByDay([outmat.accuracy outmat.actionRate outmat.stimulus outmat.bias outmat.missBias],...
    outmat.probeData,outmat.probeTrialNum,outmat.dayLen,...
    {'Accuracy','ActionRate','StimulusFreq','ActionBias','MissBias'},mouse);
    
saveas(f,[figPath sep mouse '.png']);
saveas(f,[figPath sep mouse '.m']);
close all;
end


%----------------------Code for Plotting------------------------------------------
function f = plotAnimalByDay(mat,probeData,probeTrialNum,dayLen,ylabels,mouse)

if nargin == 2; ylabels = {[],[]};end
cumsumDayLen = [0 cumsum(dayLen)];

f = figure; set(f,'Units','Normalized','OuterPosition',[0.05,0.58,0.32,0.40])
for i = 1:size(mat,2)
    subplot_tight(size(mat,2),1,i,[0.05 0.05]); hold on;
    if i == 1 || i==2; ylimm = [0 1]; else; ylimm = [-1 1]; end
    if i == 1; plot([1 cumsumDayLen(end)],[0.5 0.5],...
            'Color',[0.8 0.8 0.8],'LineWidth',2);
    elseif i==3 || i==4 || i==5; plot([1 cumsumDayLen(end)],[0 0],...
            'Color',[0.8 0.8 0.8],'LineWidth',2); 
    end
    
    for j = 1:length(dayLen)
        plot([cumsumDayLen(j+1) cumsumDayLen(j+1)], ylimm, 'Color', [0.8 0.8 0.8], 'LineWidth',1.5);
        plot(cumsumDayLen(j)+1:cumsumDayLen(j+1),mat(cumsumDayLen(j)+1:cumsumDayLen(j+1),i),...
            'Color',matlabColors(1,0.8),'LineWidth',2);             
    end   

    try
        probeData(isnan(probeData(:,1)),:) = [];
        if i == 1; plot(probeTrialNum,mean(probeData(:,1:2),2),'--o',...
                'Color',[0.2 0.2 0.2],'LineWidth',2);
        elseif i == 2; plot(probeTrialNum,mean(probeData(:,3:4),2),'--o',...
                'Color',[0.2 0.2 0.2],'LineWidth',2);
        elseif i == 4; plot(probeTrialNum,probeData(:,5),'--o',...
                'Color',[0.2 0.2 0.2],'LineWidth',2);
        end    
    catch
        disp('probe bugged')
    end

    
    %xlim([1 2400])
    xlim([1 cumsumDayLen(end)])
    ylim(ylimm);
    ylabel(ylabels{i})
    if i~= size(mat,2); xticks([]); else; xlabel('Trials'); end
    if i==1; title(mouse); end
end

end