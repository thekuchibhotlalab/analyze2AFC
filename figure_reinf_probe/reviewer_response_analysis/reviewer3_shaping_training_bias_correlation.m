%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT_bin40(mice);
mouseMega = wheel2AFCmega(allMouse);

%% 
preTrainingMove = {[52 44; 133 216], [167 230; 258 172], [123,236; 189,328], [134,208; 208,243; 233,204; 240 315; 218 244],... %062-066
    [58,84; 211,285; 166,234; 264,442; 203,319], [219,193; 237 262; 164 218; 186 272; 139 203],... %067-068
    [213 259; 282 311; 210 246; 252 335; 268 268], [130 141; 178 210; 173 235; 135 155; 150 200],... %069-107
    [114 103; 187 178; 199 256; 189 228; 243 225], [239 199; 247 238; 237 166; 294 248; 307 189],...%109-111
    [144 184; 106 109; 187 206; 194 240; 212 245], [143 164; 135 201; 136 138; 148 193; 156 220],...
    [170 165; 171 218; 173 198; 266 301; 175 214]};
trainingBiasMove = [172 450; 249 921; 281 1053; 832 766; 946 1771; 330 768; 1044 1056; 694 1657; 1047 1250;...
    2716 2083; 1311 1362; 958 873; 968 1520];
preTrainingBias = fn_cell2mat(cellfun(@(x)(mean(x(:,:)./...
    repmat(sum(x(:,:),2),[1 size(x,2)]),1)),preTrainingMove,'UniformOutput',false),1);

trainingBiasMove = trainingBiasMove ./ repmat(sum(trainingBiasMove,2),[1 2]);
preTrainingBias = preTrainingBias(:,1) - preTrainingBias(:,2);
trainingBias = [];
trainingBiasMoveBias = trainingBiasMove(:,1) - trainingBiasMove(:,2);
for i = 1:nMouse
    temp = []; 
    for j = 1:5
        tempFlag = allMouse{i}.behav.day==j; 
        [bias, acc] = fn_getAccBias(allMouse{i}.behav.stimulus(tempFlag),...
            allMouse{i}.behav.stimulus(tempFlag)== allMouse{i}.behav.action(tempFlag));
    temp(j) = bias;
    end
    trainingBias(i) = mean(temp);
end

figure; lm = fn_plotScatterCorr(preTrainingBias',trainingBias)
xlim([-0.5 0.5]); ylim([-0.5 0.5])
figure; 
fn_plotComparison({preTrainingBias,trainingBias},'paired','true');
figure; 
fn_plotComparison({abs(preTrainingBias),abs(trainingBias)},'paired','true')
xlim([0 3]); xticks([1 2]); ylim([0 0.5])
[r,p] = corr(preTrainingBias,trainingBias')
[r,p] = corr(trainingBias',trainingBiasMoveBias)
figure; scatter(trainingBias',trainingBiasMoveBias)
%% stability of wheel training

preTrainingBiasAllday = cellfun(@(x)(x ./ repmat(sum(x,2),[1 2])),preTrainingMove,'UniformOutput',false);
preTrainingBiasAllday = cellfun(@(x)(x(:,1)-x(:,2)),preTrainingBiasAllday,'UniformOutput',false);
preTrainingBiasAllday = fn_cell2matFillNan(preTrainingBiasAllday);

figure; fn_plotMeanSampleLine(1:5, abs(preTrainingBiasAllday'),...
    {'LineWidth',2,'Color',[0.2 0.2 0.2]},{'Color',[0.8 0.8 0.8]});
ylim([0 0.4])

preTrainingBiasAllday10 = preTrainingBiasAllday(:,4:end);
preTrainingBiasAllday10(preTrainingBiasAllday10>0.15) = 1; preTrainingBiasAllday10(preTrainingBiasAllday10<-0.15) = -1;
preTrainingBiasAllday10(preTrainingBiasAllday10~=1 & preTrainingBiasAllday10~=-1) = 0;
preTrainingBiasAllday10Trans = sum(diff(preTrainingBiasAllday10)~=0);

figure; histogram(preTrainingBiasAllday10Trans,-0.5:1:5.5);
xlim([-1 6]); ylim([0 10])