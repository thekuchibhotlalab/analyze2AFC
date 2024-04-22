%% deterministic choice

clear;
startingFactor = 1; nTrials = 20; updateFactor= 1; forgetFactor = 1;
param{1} = initializeParam(startingFactor,nTrials); % param for L stim
param{2} = initializeParam(startingFactor,nTrials); % param for R stim

choices = nan(1,nTrials); outcomPred = {}; policy = nan(1,nTrials);
stim = nan(1,nTrials);
for t = 1:nTrials
    currStim = double(rand()>0.5);stim(t) = currStim;
    [ tempChoice, outcome,outcomPred{t},policy(t)] = generateChoice(param,currStim,t,updateFactor,forgetFactor,'');

    param = update(param,currStim, tempChoice, outcome,t,updateFactor,forgetFactor);
    choices(t) = tempChoice; 
end

%% figures

figure; subplot(1,2,1); plot(param{1}.mu); ylim([0 1])
subplot(1,2,2); plot(param{2}.mu); ylim([0 1])

figure; subplot(1,2,1); plot(param{1}.epi); ylim([0 1])
subplot(1,2,2); plot(param{2}.epi); ylim([0 1])

figure; plot(policy)

figure; subplot(1,3,1); plot(outcomPred{2});
subplot(1,3,2); plot(outcomPred{20})
subplot(1,3,3); plot(outcomPred{10})
%% probabilisti choice
clear;
startingFactor = 1; nTrials = 100; updateFactor= 1; forgetFactor = 1;
param{1} = initializeParam(startingFactor,nTrials); % param for L stim
param{2} = initializeParam(startingFactor,nTrials); % param for R stim

choices = nan(1,nTrials); outcomPred = {}; policy = nan(1,nTrials);
stim = nan(1,nTrials);
for t = 1:nTrials
    currStim = double(rand()>0.5);stim(t) = currStim;
    [ tempChoice, outcome,outcomPred{t},policy(t)] = generateChoice(param,currStim,t,updateFactor,forgetFactor,'prob');

    param = update(param,currStim, tempChoice, outcome,t,updateFactor,forgetFactor);
    choices(t) = tempChoice; 
end
%% figures
figure; subplot(1,2,1); plot(param{1}.mu); ylim([0 1])
subplot(1,2,2); plot(param{2}.mu); ylim([0 1])

figure; subplot(1,2,1); plot(param{1}.epi); ylim([0 1])
subplot(1,2,2); plot(param{2}.epi); ylim([0 1])


figure; plot(policy)

figure; subplot(1,3,1); plot(outcomPred{2});
subplot(1,3,2); plot(outcomPred{20})
subplot(1,3,3); plot(outcomPred{100})
%% all functions

function param = initializeParam(startingFactor,nTrials)
param.alpha = nan(nTrials+1,2); param.alpha(1,:) = [1 1] * startingFactor; 
param.beta = nan(nTrials+1,2); param.beta(1,:) = [1 1] * startingFactor; 
param.mu = nan(nTrials+1,2); param.mu(1,:) = param.alpha(1,:) ./ (param.alpha(1,:) + param.beta(1,:));
param.var = nan(nTrials+1,2); param.var(1,:) = param.alpha(1,:) .* param.beta(1,:) ./ (( param.alpha(1,:)+param.beta(1,:) +1) .* ( param.alpha(1,:)+param.beta(1,:) ).^2);
param.epi = nan(nTrials+1,2); param.epi(1,:) = param.var(1,:)./ (param.mu(1,:).*(1-param.mu(1,:)));
end

function [choice, outcome,outcomPred,bestPolicy] = generateChoice(param,currStim, t,updateFactor,forgetFactor,probFlag)
horizon = 5; 
actionSeq = zeros(2^horizon,horizon);  actionSeq = initSeq(actionSeq);
stimRewSeq = zeros(2^(2*horizon-1),2*horizon-1); stimRewSeq = initSeq(stimRewSeq);
stimRewSeq = cat(2, ones(2^(2*horizon-1),1)*currStim, stimRewSeq);
%probSeq = zeros(size(stimRewSeq,1),2^horizon);

tempStateProb = nan(size(stimRewSeq,1),size(actionSeq,1),horizon);
tempRewHorizon = nan(size(stimRewSeq,1),size(actionSeq,1),horizon);
for j = 1:(size(stimRewSeq,1))
    tempStimRew = stimRewSeq(j,:);
    tempStim = tempStimRew(1:length(tempStimRew)/2); 
    tempRew = tempStimRew(length(tempStimRew)/2+1:end);
    tempRewHorizon(j,:,:) = repmat(tempRew,[size(actionSeq,1) 1]); 
    for i = 1:size(actionSeq,1)
        %tempRewardProb = nan(1,horizon);
        tempParam = param;
        for l = 1:horizon
            tempRewardProb = tempParam{tempStim(l)+1}.mu(t+l-1,actionSeq(i,l)+1);
            tempStateProb(j,i,l) = (0.5 + double(l==1)*0.5) * (tempRew(l)*tempRewardProb + (1-tempRew(l))*(1-tempRewardProb));
            %tempParam = update(tempParam,tempStim(l),actionSeq(i,l),tempRew(l),t+l-1,updateFactor,forgetFactor);
            %tempRewardProb(l) = tempParam{tempStim(l)+1}.mu(t+l-1,actionSeq(i,l)+1);
            tempParam = update(tempParam,tempStim(l),actionSeq(i,l),tempRew(l),t+l-1,updateFactor,forgetFactor);
            
        end
        %probSeq(j,i) = sum(tempRewardProb .* [1 1 1 1 1]); 
    end
end
probSeq = prod(tempStateProb,3);
weights = [1;1;1;1;1]; rewHorizon = reshape(tempRewHorizon, [size(stimRewSeq,1)*size(actionSeq,1),horizon]);
rewHorizon = rewHorizon * weights; rewHorizon = reshape(rewHorizon, [size(stimRewSeq,1) size(actionSeq,1)]);
probSeq = probSeq .* rewHorizon;
outcomPred = sum(probSeq,1);

if strcmp(probFlag,'prob')
 bestPolicy = randsample(length(outcomPred),1,true,softmax(outcomPred'));
else
 [~,bestPolicy] = max(outcomPred);
end
%if t<20; choice = 0; else; end
choice = actionSeq(bestPolicy,1);
outcome = double(currStim==choice);
if t==10
    disp('10')
end



function seq = initSeq(seq)
    for k = 0:(size(seq,1)-1)
        a = cellfun(@str2num,num2cell(dec2bin(k)));
        seq(k+1,size(seq,2)-length(a)+1:end) = a;
    end
end

end

function param = update(param, stim, action, outcome,t,updateFactor,forgetFactor)
    % stim presented, update
    param{stim+1}.alpha(t+1,:)= [param{stim+1}.alpha(t,1)* forgetFactor+ (1-action)*outcome*updateFactor, ...
        param{stim+1}.alpha(t,2)* forgetFactor+ action*outcome*updateFactor];
    param{stim+1}.beta(t+1,:) = [param{stim+1}.beta(t,1)* forgetFactor + (1-action)*(1-outcome)*updateFactor, ...
        param{stim+1}.beta(t,2)* forgetFactor + action*(1-outcome)*updateFactor];

    param{stim+1}.mu(t+1,:) = param{stim+1}.alpha(t+1,:) ./ (param{stim+1}.alpha(t+1,:) + param{stim+1}.beta(t+1,:));
    param{stim+1}.var(t+1,:) = param{stim+1}.alpha(t+1,:) .* param{stim+1}.beta(t+1,:) ./ ((param{stim+1}.alpha(t+1,:) + param{stim+1}.beta(t+1,:)+1) .*...
        ( param{stim+1}.alpha(t+1,:)+ param{stim+1}.beta(t+1,:) ).^2);
    param{stim+1}.epi(t+1,:) = param{stim+1}.var(t+1,:)./ (param{stim+1}.mu(t+1,:).*(1-param{stim+1}.mu(t+1,:)));

    % stim not presented, no update
    param{2-stim}.alpha(t+1,:) = param{2-stim}.alpha(t,:);
    param{2-stim}.beta(t+1,:) = param{2-stim}.beta(t,:);
    param{2-stim}.mu(t+1,:) = param{2-stim}.mu(t,:);
    param{2-stim}.var(t+1,:) = param{2-stim}.var(t,:);
    param{2-stim}.epi(t+1,:) = param{2-stim}.epi(t,:);
end