%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT(mice);
mouseMega = wheel2AFCmega(allMouse);
%%
reinfThre = 0.60:0.02:0.8;
avgStd = nan(size(reinfThre));
for i = 1:length(reinfThre)
    try
        [avgStd(i)] = alignStuff(reinfThre(i),mouseMega);
    catch
    end
end
%%
figure; plot(reinfThre,avgStd,'Color',[0.2 0.2 0.2]); hold on; scatter(reinfThre,avgStd,30,[0.2 0.2 0.2],'filled');
yticks(0.072:0.004:0.08); ylim([0.072 0.08])
%%
function [avgStd] = alignStuff(reinfThre,mouseMega)


% align data by minimizing performance variance
probeThre = nan; probeTrialBin = 400;
outCell = mouseMega.objFun('binProbeByTrialFromLearningOnset',{[reinfThre probeThre],probeTrialBin});


reinfData = mouseMega.loadReinf;
reinfAlignPoint = cell2mat(outCell{3});
[reinfDataAlign,reinfAlignPoint] = attachNan(reinfData, reinfAlignPoint);


a = sum(isnan(reinfDataAlign.acc),2); tempStart = find(a==5,1,'first'); tempEnd = find(a==5,1,'last');

avgStd = mean(nanstd(reinfDataAlign.acc(tempStart:tempEnd,:),0,2));
end

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