function [wheelDownSample,badFlag] = getDownSample(mouseMega,idx)
if isempty(idx); idx = 1:mouseMega.nMouse; end
wheelDownSample = cell(1,length(idx));
for i = 1:length(idx)
    [wheelPos, badFlag{i}] = fn_removeBadTrials(mouseMega.mouseCell{idx(i)}.behavVar.wheelPos_aligned);
    
    % normalize the max position of the wheel to 1; downsample every 50ms 
    downSampleRate = 50; tempTimeDownSampleIdx = 2001:downSampleRate:4501;
    tempMax = wheelPos(:,tempTimeDownSampleIdx(1):tempTimeDownSampleIdx(end)); tempMax = max(tempMax(:));
    wheelPos = wheelPos /tempMax;
    tempWheelDownSample = wheelPos(:,tempTimeDownSampleIdx);
    wheelDownSample{i} = tempWheelDownSample;
end 


end