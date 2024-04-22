clear;
%mice = {'zz054','zz063','zz066','zz067','zz068','zz069'};
%mice = {'zz066','zz067','zz068','zz069'};
%mice = {'zz081','zz082','zz083'}; %FM mice FM_UpDown
%mice = {'zz075','zz076','zz077','zz048','zz049'}; %FM 2 oct/s FM_Multi_Oct
%mice = {'zz070','zz071','zz072','zz073'}; %FM_Oct FM_Oct_prob
%mice = {'zz065','zz066','zz079'}; %AM_puretone
%mice = {'zz060','zz062','zz063','zz064'}; %AM_WN
%mice = {'zz101','zz102','zz103','zz104',...
%    'zz105','zz106','zz107','zz108','zz109','zz110','zz111','zz112',...
%    'zz113','zz114','zz115','zz116','zz117','zz118','zz119','zz120',...
%    'zz121','zz122','zz123','zz124','zz125'};
%mice = {'msb01','msb02','msb03','msb04'};
%mice = {'msb08','msb09'};

%mice = {'zz107','zz109','zz111','zz112','zz113','zz115'};
%mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
mice = {'sz06','sz07','sz08','sz09','sz10'};
%% save coulbourn data to datapath
cellfun(@readCoulbournData,mice,'UniformOutput',false);
%%
error('stop!')
%% plot learning curve of each mouse
cellfun(@(x)plotLearningCurveTrial(x,{'puretone'}),mice,'UniformOutput',false);
%% align mouse learning curve with miss trials
plotLearningCurveTrialAll(mice,{'puretone'});
%% align mouse learning curve without miss trials
plotLearningCurveTrialAllNoMiss(mice,{'puretone'});
%plotLearningCurveTrialAllNoMiss(mice,{'FM_UpDown'});
%plotLearningCurveTrialAllNoMiss(mice,{'FM_Multi_Oct','FM_Multi_Oct_Prob'});
%plotLearningCurveTrialAllNoMiss(mice,{'AM_puretone'});
%plotLearningCurveTrialAllNoMiss(mice,{'AM'});
%% PLOT THIS COHORT
mice = {'zz111','zz115'};
cellfun(@(x)plotLearningCurveTrial(x,{'FM_Dir_pip'}),mice,'UniformOutput',false);
%% PLOT THIS COHORT
mice = {'zz121','zz122','zz123','zz124','zz125'};
cellfun(@(x)plotLearningCurveTrial(x,{'FM_Dir_pip'}),mice,'UniformOutput',false);

%% PLOT THIS COHORT
clear;
mice = {'zz116','zz117','zz118','zz119','zz120'};
cellfun(@(x)plotLearningCurveTrial(x,{'PT_pip'}),mice,'UniformOutput',false);
clear;
mice = {'zz116','zz117','zz118','zz119','zz120'};
cellfun(@(x)plotLearningCurveTrial(x,{'FM_Dir_pip'}),mice,'UniformOutput',false);
%% align mouse learning curve without miss trials
[tempSum,blockLen,midTransFreq] = plotBiasTransition(mice,{'puretone'});

%%
clear;
mice = {'zz080PT','zz081PT','zz082PT','zz083PT','zz100','zz097','zz098','zz099','zz101','zz102','zz103','zz104','zz105'...
    'zz106','zz107','zz108','zz109','zz110'};
cellfun(@readCoulbournData,mice,'UniformOutput',false);
%%
clear;
mice = {'zz101','zz102','zz104','zz105'};
cellfun(@(x)plotLearningCurveTrial(x,{'FM_Dir_Dur_pip'}),mice,'UniformOutput',false);
%%
clear;
mice = {'zz100','zz097','zz098','zz099'};
cellfun(@(x)plotLearningCurveTrial(x,{'PT_pip'}),mice,'UniformOutput',false);
mice = {'zz097','zz098','zz099','zz100'};
cellfun(@(x)plotLearningCurveTrial(x,{'PT_pip'}),mice,'UniformOutput',false);
mice = {'zz101','zz102','zz103','zz104','zz105'};
cellfun(@(x)plotLearningCurveTrial(x,{'FM_Dir_Dur_pip'}),mice,'UniformOutput',false);
mice = {'zz106','zz107','zz108','zz109','zz110'};
cellfun(@(x)plotLearningCurveTrial(x,{'puretone'}),mice,'UniformOutput',false);

%% align mouse learning curve with miss trials
clear; 
mice = {'zz097','zz098'};
%plotLearningCurveTrialAllNoMiss(mice,{'PT_pip'});
plotBiasTransition(mice,{'PT_pip'});

%% align mouse learning curve with miss trials
clear; 
mice = {'zz101','zz102','zz105'};
%plotLearningCurveTrialAllNoMiss(mice,{'FM_Dir_Dur_pip'});
plotBiasTransition(mice,{'FM_Dir_Dur_pip'});
%% align mouse learning curve with miss trials
clear; 
mice = {'zz107','zz109'};
%plotLearningCurveTrialAllNoMiss(mice,{'puretone'});
plotBiasTransition(mice,{'puretone'});

%% align mouse learning curve with miss trials
clear; 
mice = {'zz107'};
plotLearningCurveTrial('zz107',{'puretone'});
%plotBiasTransition(mice,{'puretone'});
%%
clear;
mice = {'zz054','zz063','zz066','zz067','zz068','zz069'};
[tempSum1,blockLen1,midTransFreq1,bins1] = plotBiasTransition(mice,{'puretone'});
close all;
mice = {'zz097','zz098'};
[tempSum2,blockLen2,midTransFreq2,bins2] = plotBiasTransition(mice,{'PT_pip'});
close all;
mice = {'zz107','zz109'};

[tempSum3,blockLen3,midTransFreq3,bins3] =plotBiasTransition(mice,{'puretone'});
close all;
mice = {'zz101','zz102','zz105'};
[tempSum4,blockLen4,midTransFreq4,bins4] =plotBiasTransition(mice,{'FM_Dir_Dur_pip'});
close all;

%%
bins = bins1;
tempSum = cat(1,tempSum1,tempSum2,tempSum3,tempSum4);
blockLen = cat(1,blockLen1,blockLen2,blockLen3,blockLen4);
midTrans = cat(1,midTransFreq1,midTransFreq2,midTransFreq3,midTransFreq4);
figure; subplot(2,1,1); hold on;

tempSum = tempSum ./ mean(mean(tempSum(:,1:2))); blockLen = blockLen ./ mean(mean(blockLen(:,1:2)));
scatter(bins(1:end-1)+bins(2)/2,nanmean(tempSum,1),20,matlabColors(1),'filled')

scatter(bins(1:end-1)+bins(2)/2,nanmean(blockLen,1),20,matlabColors(6),'filled')
f = lsline;
set(f(1), 'Color', matlabColors(6),'LineWidth',2); set(f(2), 'Color', matlabColors(1),'LineWidth',2)
legend(f,{'Bias Block Len','Bias Block Total'})
ylabel('Normalized Unit'); xlabel('Trials');ylim([0.2 1.4]); yticks([0.2:0.4:1.4]); xlim([0 2500])

subplot(2,1,2); hold on;
midTransFreq = histcounts(midTrans,bins)./sum(~isnan(tempSum),1);
scatter(bins(1:end-1)+bins(2)/2,midTransFreq./mean(midTransFreq(1:2)),20,matlabColors(2),'filled');
f = lsline;

set(f(1), 'Color', matlabColors(2),'LineWidth',2)
legend(f,{'Transitions'})
ylabel('Normalized Unit'); xlabel('Trials');ylim([0.5 1.8]); yticks([0.2:0.4:1.8]); xlim([0 2500])
