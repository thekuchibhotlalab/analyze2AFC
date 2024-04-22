%% get data
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
allMouse = fn_getObjPT_bin40(mice);
mouseMega = wheel2AFCmega(allMouse);

%% get the trial number of transition each day 
biasFirstHalf = {}; biasLastHalf = {};
for i = 1:length(allMouse)
    tempBiasFirst = []; tempBiasLast = [];
    for j = 1:max(allMouse{i}.behav.day)
        tempDayIdx = allMouse{i}.behav.day==j;
        tempDayStart = find(tempDayIdx,1,'first');
        tempDayEnd = find(tempDayIdx,1,'last');
        tempDayMid = tempDayStart + round((tempDayEnd-tempDayStart+1)/2)-1;

        tempBiasFirst(j) = nanmean(abs(allMouse{i}.behav.bias(tempDayStart+20:tempDayMid-20)));
        tempBiasLast(j) = nanmean(abs(allMouse{i}.behav.bias(tempDayMid+20:tempDayEnd-20)));

    end
    biasFirstHalf{i} = tempBiasFirst; 
    biasLastHalf{i} = tempBiasLast; 
end

biasFirstHalf = fn_cell2matFillNan(biasFirstHalf);
biasLastHalf = fn_cell2matFillNan(biasLastHalf);

%% make plot across days
figure; hold on;
tempFirst = nanmean(biasFirstHalf(:,1:10),1);
tempLast = nanmean(biasLastHalf(:,1:10),1);
fn_plotMeanErrorbar(1:size(biasFirstHalf,2),biasFirstHalf,matlabColors(1),matlabColors(1),...
    {}, {'faceAlpha',0.2})

fn_plotMeanErrorbar(1:size(biasLastHalf,2),biasLastHalf,matlabColors(2),matlabColors(2),...
    {}, {'faceAlpha',0.2})
xlim([1 10])
figure; 
fn_plotComparison({tempFirst,tempLast},'paired',true);

p = [];
for i = 1:10
    [p(i)] = signrank(biasFirstHalf(:,i),biasLastHalf(:,i));

end
%% make plot by animal
tempFirst = nanmean(biasFirstHalf(:,1:10),2);
tempLast = nanmean(biasLastHalf(:,1:10),2);
figure; hold on;
fn_plotComparison({tempFirst,tempLast},'paired',true)