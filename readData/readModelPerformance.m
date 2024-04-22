function weightCell =readModelPerformance(mice)
dataPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData';
sep = '\';
plotLim = [1 2400];
%mice = {'zz054','zz062','zz063','zz066','zz067','zz068','zz069'};

%% 1.1 - MAKE ACCURACY LEARNING CURVE ACROSS MICE
learningCurveBin = 50;
accuracyCell = {};
for i = 1:length(mice)
    mouse = mice{i};
    if exist([dataPath sep mouse '.mat'])
        behavior = load([dataPath sep mouse '.mat']);
        accuracyCell{i} = smoothdata(behavior.correct,'movmean', learningCurveBin);
    end
end

if isempty(accuracyCell); weightCell = {}; return; end

%% 1.2 - MAKE REGRESSION LEARNING CURVE ACROSS MICE
weightCell = {};
for i = 1:length(mice)
    mouse = mice{i};
    model = load([dataPath sep mouse 'psytrack_SBAArArs.mat']);
    maxAcc = max(accuracyCell{i});
    minAcc = min(accuracyCell{i});
    weightCell{i} = model.wMode(5,:);
    weightCell{i} = (1+exp(-weightCell{i})).^(-1);
    %weightCell{i} = (weightCell{i}-min(weightCell{i})) / (max(weightCell{i})...
    %    - min(weightCell{i})) * (maxAcc-minAcc)+minAcc;
end

end