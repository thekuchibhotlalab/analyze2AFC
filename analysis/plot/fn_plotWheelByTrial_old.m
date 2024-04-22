function fn_plotWheelByTrial(mouseMega,exampleAnimal,trialnum,varargin)

p = inputParser;
p.KeepUnmatched = true;
% ------ Varargin for superclass
p.addParameter('probeIdx', []);
p.parse(varargin{:});

tempWheel = mouseMega.mouseCell{exampleAnimal}.behavVar.wheelPos_aligned(trialnum,:);
tempResp = mouseMega.mouseCell{exampleAnimal}.behav.responseType(trialnum);
tempStim = mouseMega.mouseCell{exampleAnimal}.behav.stimulus(trialnum);
tempRT = mouseMega.mouseCell{exampleAnimal}.behav.reactionTime(trialnum);
downSample = 10; tempRT = round(tempRT*1000/downSample);
tempSpeed = diff(tempWheel,1,2); tempSpeed = fn_binSum(tempSpeed, downSample,2);
nSample = round(tempSpeed/downSample);

hold on;
idxC1 = {};idxC2 = {}; idxIC1 = {};idxIC2 = {};
for i = 1:size(tempSpeed,1) % loop through all 
if tempStim(i) == 1
    [~,tempC]= find(tempSpeed(i,:)>1); idxC1{i} = ones(1,length(tempC))*i; idxC2{i} = tempC;
    [~,tempIC]= find(tempSpeed(i,:)<-1); idxIC1{i} = ones(1,length(tempIC))*i; idxIC2{i} = tempIC;
    
else 
    [~,tempC]= find(tempSpeed(i,:)<-1); idxC1{i} = ones(1,length(tempC))*i; idxC2{i} = tempC;
    [~,tempIC]= find(tempSpeed(i,:)>1); idxIC1{i} = ones(1,length(tempIC))*i; idxIC2{i} = tempIC;
end 
end

idxC1 = fn_cell2mat(idxC1,2);idxC2 = fn_cell2mat(idxC2,2);
idxIC1 = fn_cell2mat(idxIC1,2);idxIC2 = fn_cell2mat(idxIC2,2);

tempArea = [];
for i = 1:length(tempRT); tempArea = cat(2, tempArea, [ cat(2,tempRT(i)+200,tempRT(i)+200);cat(2,i-1,i)]); end
tempArea = cat(2,tempArea,cat(1,200*ones(1,length(tempRT)+1),length(tempRT):-1:0));
fill(tempArea(1,:),tempArea(2,:),[0.8 0.8 0.8],'LineStyle','None');

if ~isempty(p.Results.probeIdx)
    probeIdx = [p.Results.probeIdx(1)-1 p.Results.probeIdx];
    probeArea = tempArea; probeFlag = probeArea(2,:);
    probeFlag = fn_findAny(probeFlag,probeIdx);
    probeArea = probeArea(:,probeFlag);
    fill(probeArea(1,:),probeArea(2,:),fn_wheelColorsPT('Probe',0.3),'LineStyle','None');
end

scatter(idxC2,idxC1-0.5,5,fn_wheelColorsPT('correct'),'filled')
scatter(idxIC2,idxIC1-0.5,5,fn_wheelColorsPT('incorrect'),'filled')
xlim([150 450]); xticks(150:50:450); xticklabels(-0.5:0.5:2.5); ylim([0 length(tempRT)])

end