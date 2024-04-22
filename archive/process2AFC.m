clear;
mice = {'zz054','zz063','zz066','zz067','zz068','zz069'};
%mice = {'zz066','zz067','zz068','zz069'};
%mice = {'zz081','zz082','zz083'}; %FM mice FM_UpDown
%mice = {'zz075','zz076','zz077','zz048','zz049'}; %FM 2 oct/s FM_Multi_Oct
%mice = {'zz070','zz071','zz072','zz073'}; %FM_Oct FM_Oct_prob
%mice = {'zz065','zz066','zz079'}; %AM_puretone
%mice = {'zz060','zz062','zz063','zz064'}; %AM_WN
%% save coulbourn data to datapath
cellfun(@readCoulbournData,mice,'UniformOutput',false);

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

%% align mouse learning curve without miss trials
plotBiasTransition(mice,{'puretone'});

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
plotLearningCurveTrialAllNoMiss(mice,{'PT_pip'});
plotBiasTransition(mice,{'PT_pip'});

%% align mouse learning curve with miss trials
clear; 
mice = {'zz101','zz105'};
plotLearningCurveTrialAllNoMiss(mice,{'FM_Dir_Dur_pip'});
plotBiasTransition(mice,{'FM_Dir_Dur_pip'});
%% align mouse learning curve with miss trials
clear; 
mice = {'zz107','zz109'};
plotLearningCurveTrialAllNoMiss(mice,{'puretone'});
plotBiasTransition(mice,{'puretone'});

%% align mouse learning curve with miss trials
clear; 
mice = 'zz097';
plotLearningCurveTrial(mice,{'FM_Dir_pip'});