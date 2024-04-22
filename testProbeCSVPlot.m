clear;
sep = '\';
loadPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\probeCSV';
mouse = 'zz076';
csvData = readtable([loadPath sep mouse '.csv']);

day = csvData.Day;

for i = 1:length(day)
    tempL = (csvData.PLCorr + csvData.PLInCorr + csvData.PLInMiss);
    tempR = (csvData.PRCorr + csvData.PRInCorr + csvData.PRInMiss);
    if  tempL> 10; csvData.PLCorr = csvData.PLCorr - tempL + 10; end
    if  tempR> 10; csvData.PRCorr = csvData.PRCorr - tempR + 10; end
    
end
accL = (csvData.PLCorr)./(csvData.PLCorr + csvData.PLInCorr);
accR = (csvData.PRCorr)./(csvData.PRCorr + csvData.PRInCorr);

missL = (csvData.PLCorr + csvData.PLInCorr)./(csvData.PLCorr + csvData.PLInCorr + csvData.PLInMiss);
missR = (csvData.PRCorr + csvData.PRInCorr)./(csvData.PRCorr + csvData.PRInCorr + csvData.PRInMiss);


nanFlag = isnan(sum([accL accR],2));
fn_figureSmartDim('hSize',0.23,'widthHeightRatio',1.5); hold on;
%plot(day(~nanFlag),accL(~nanFlag),'-o','LineWidth',2);
plot(day(~nanFlag),accR(~nanFlag),'-o','LineWidth',2,'Color',matlabColors(2));
plot([day(1) day(end)],[0.5 0.5],'Color',[0.8 0.8 0.8],'LineWidth',2);
xlabel('Day'); ylabel('Accuracy'); xlim([day(1) day(end)])
fn_figureSmartDim('hSize',0.23,'widthHeightRatio',1.5); hold on;
%plot(day(~nanFlag),missL(~nanFlag),'-o','LineWidth',2);
plot(day(~nanFlag),missR(~nanFlag),'-o','LineWidth',2,'Color',matlabColors(2));
xlabel('Day'); ylabel('Action Rate'); xlim([day(1) day(end)])

%%
mice = { 'zz075', 'zz076', 'zz077', 'zz081', 'zz082', 'zz083'};
[accL_All, accR_All,missL_All,missR_All] = cellfun(@testProbeCSV,mice,'UniformOutput',false);
accL_All = fn_cell2mat(accL_All,1); accR_All = fn_cell2mat(accR_All,1);
missL_All = fn_cell2mat(missL_All,1); missR_All = fn_cell2mat(missR_All,1);

figure; errorbar(nanmean(accL_All,1), nanstd(accL_All,0,1)./sqrt(sum(~isnan(accL_All),1)),'LineWidth',2); 
hold on; errorbar(nanmean(accR_All,1), nanstd(accR_All,0,1)./sqrt(sum(~isnan(accR_All),1)),'LineWidth',2); 
xlim([0.5 3.5]);ylim([0.3 1]);ylabel('Accuracy'); xticks([1 2 3]); xticklabels({'Pre','Acquisition','Post'})

figure; errorbar(nanmean(missL_All,1), nanstd(missL_All,0,1)./sqrt(sum(~isnan(missL_All),1)),'LineWidth',2); 
hold on; errorbar(nanmean(missR_All,1), nanstd(missR_All,0,1)./sqrt(sum(~isnan(missR_All),1)),'LineWidth',2); 
xlim([0.5 3.5]);ylim([0.45 1]);ylabel('Action Rate'); xticks([1 2 3]); xticklabels({'Pre','Acquisition','Post'})
