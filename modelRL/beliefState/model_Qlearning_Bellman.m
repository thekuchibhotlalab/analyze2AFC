%% initiate model (i.e. initiate alpha and beta to 0) and stimuli train for all trials
nTrials = 20; n = 2;
model = initiateModel();
stimuli = generateStimuliTrain(nTrials);
modelList = {}; choice = []; Q = [];
for trial = 1:nTrials
    % For each trial, to make a choice, we fist evaluate horizon (here use ‘n’ to evaluate the outcome of future n action, thus (n-1) step horizon beyond the current step)
    [maxAction,Q(trial)] = evaluateHorizon(model, n);
    % maxAction is a 1-by-2 array of the best action for each current state. Select the best action based on state
    tempChoice = maxAction(stimuli(trial)+1);
    % compute reward
    if stimuli(trial) == tempChoice; reward=1; else; reward = 0; end 
    % update model
    [model] = updateBelief(model, stimuli(trial), tempChoice,reward); 
    modelList{trial} = model; 
    choice(trial) = tempChoice;
end 

figure; plot(Q); xlabel('trials'); ylabel('Q')


%% function to evaluate horizon
function [maxAction,bestQ] = evaluateHorizon(model, n)
tau = 1;
% Generate all stimulus sequences, action sequences, and the probability of getting reward on each time on each sequence. 1st dim: all possible stim/act/reward sequences; 2nd dimension: identity of each event sequence binarized as 0/1, with all time concatenated together. E.g. [0 1 0, 1 0 0] means: at time = t, stim=0, action=1, reward=0.
allSeq = zeros (2^(3*n), 3*n);
allSeq = initSeq(allSeq);tempSeq = allSeq;
% Loop through (n-1)-step horizon, total of n timestep including the current time step. Start optimizing from last time step as in Bellman equation
for time = n:-1:1 
    % we will need to marginalize/sum up over all possible states at time=n, so we first find all the sequences (by finding the ‘unique sequences’ at time 1 to n-1) 
    [margSeq, margSeqIdx]= unique(tempSeq(:,1:3*(time-1)), 'rows');
     bestQ = []; 
    if ~isempty(margSeq) % earlier time steps
        for seq = 1:size(margSeq,1)
            % update the belief given each past (state,action,reward) sequence
            newBeliefModel = updateBelief(model,margSeq(seq,1:3:end),margSeq(seq,2:3:end),margSeq(seq,3:3:end));
    
            % evaluate immediate reward gain at current time step. First: get the current (stim,action,reward) sequence. This will be all the possible sequences that comes after the sequence that we need to marginalize
            tempIdx = sum(tempSeq(:,1:3*(time-1)) == repmat(margSeq(seq,:),[size(tempSeq,1) 1]),2)==3;
            currSeq = tempSeq(tempIdx,3*(time-1)+1: 3*(time));
            % evaluate the gain of each sequence according to belief model.
            currRew = weighCurrRew(newBeliefModel,currSeq); currRew(isnan(currRew)) = 0; currRew(isinf(currRew)) = 0; 
            % evaluate future reward, if we are not at the last step of horizon. future is the value passed from the previous iteration, representing the maximizing Q value of the next time step, summed across state types
            if time~=n; futureReward = weighFutureRew(newBeliefModel, currSeq, futureQ(tempIdx)); 
            else; futureReward = 0; 
            end
            % sum up reward. Now, currQ should be a list of Q value of each (s,a) pair at current time step
            currQ = currRew + tau*futureReward;
            % for each possible current state (i.e. stimulus 1 or stimulus 2), find action that maximize Q, and sum the best Q across states
            [maxQ, maxAction] = selectBestPolicyForEachState(currQ, currSeq);
            bestQ(seq) = maxQ;     
        end
        tempSeq = margSeq;
        % the best Q for each marginalized state sequence has been determined. This Q becomes the future Qe when we evaluate the next step
        futureQ = bestQ;
    else  % when time point reach 1 (current time point), the maximizing action for each state should become the choice of the animal (or maybe pass through a logistic choice function)
        currSeq = tempSeq;
        currRew = weighCurrRew(model,currSeq); currRew(isnan(currRew)) = 0; currRew(isinf(currRew)) = 0; 
        futureReward = weighFutureRew(model, currSeq, futureQ);
        currQ = currRew + tau*futureReward;
        [maxQ, maxAction] = selectBestPolicyForEachState(currQ, currSeq);
        bestQ = maxQ;  

    end
    
 
end
end

%% function of updating belief
function model= updateBelief(model, s,a,reward)
    nTimePoint = length(s); w = 1; 
    for time = 1:nTimePoint  %update the belief about the model starting from time 
        if reward(time)==1
        model.alpha(s(time)+1,a(time)+1) = w*model.alpha(s(time)+1,a(time)+1)+1; %w is the forgetting factor from 0 to 1
	    else % not rewarded
        model.beta(s(time)+1,a(time)+1) = w*model.beta(s(time)+1,a(time)+1)+1; 
        end
    end
end

function model= initiateModel()
    model = struct(); 
    model.alpha = ones(2,2); %w is the forgetting factor from 0 to 1
    model.beta = ones(2,2); 
end

% function of selecting best policy 
%function   [maxQ, maxAction] = selectBestPolicyForEachState(currQ, currSeq)
%maxIdx = argmax(currQ); % find the policy that give the max Q
%maxQ = currQ(maxIdx); maxAction = currSeq (maxIdx,2); % take the best Q value and the action
%end 

function stimuli = generateStimuliTrain(nTrials)
    stimuli = double(rand(1,nTrials)>0.5);

end

function seq = initSeq(seq)
    for k = 0:(size(seq,1)-1)
        a = cellfun(@str2num,num2cell(dec2bin(k)));
        seq(k+1,size(seq,2)-length(a)+1:end) = a;
    end
end

function currRew = weighCurrRew(model,currSeq)
    tempRew = zeros(size(currSeq,1),1);
    tempRew(currSeq(:,1)==0 & currSeq(:,2)==0 & currSeq(:,3)==1) = model.alpha(1,1) / (model.alpha(1,1)+model.beta(1,1));
    tempRew(currSeq(:,1)==0 & currSeq(:,2)==1 & currSeq(:,3)==1) = model.alpha(1,2) / (model.alpha(1,2)+model.beta(1,2));
    tempRew(currSeq(:,1)==1 & currSeq(:,2)==0 & currSeq(:,3)==1) = model.alpha(2,1) / (model.alpha(2,1)+model.beta(2,1));
    tempRew(currSeq(:,1)==1 & currSeq(:,2)==1 & currSeq(:,3)==1) = model.alpha(2,2) / (model.alpha(2,2)+model.beta(2,2));
    currRew = tempRew;

end


function [maxQ, maxAction] = selectBestPolicyForEachState(currQ, currSeq)
    seqL = currSeq(currSeq(:,1)==0,:);
    [maxQL, maxIdx] = max(currQ(currSeq(:,1)==0));
    maxAction(1) = seqL(maxIdx,2);
    seqR = currSeq(currSeq(:,1)==1,:);
    [maxQR, maxIdx] = max(currQ(currSeq(:,1)==1));
    maxAction(2) = seqR(maxIdx,2);
    maxQ = (maxQL+maxQR)/2;
end


function  futureReward = weighFutureRew(model, currSeq, futureQ)
    for i = 1:size(currSeq,1)
        tempProb = model.alpha(currSeq(i,1)+1,currSeq(i,2)+1)/...
            (model.alpha(currSeq(i,1)+1,currSeq(i,2)+1)+model.beta(currSeq(i,1)+1,currSeq(i,2)+1));
        stateProb(i) = currSeq(i,3) * tempProb + (1-currSeq(i,3)) * (1-tempProb);
    end
    futureReward = sum(stateProb .* futureQ) / sum(stateProb);
end
