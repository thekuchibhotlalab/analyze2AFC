clear;
startingFactor = 1; updateFactor = 0.05; 
param.alpha = [1 1] * startingFactor; param.beta = [1 1] * startingFactor; 
param.mu = param.alpha ./ (param.alpha + param.beta);
param.var = param.alpha .* param.beta ./ (( param.alpha+param.beta +1) .* ( param.alpha+param.beta ).^2);
param.epi = param.var./ (param.mu.*(1-param.mu));
nTrials = 10000;
for t = 1:nTrials
    %stim = rand()>0.5; stim = stim+1; 
    [choiceP, choice, outcome] = generateChoice(param,t);
    param = update(param, choice, outcome,t,updateFactor);

end
%% plot things
figure; 
subplot(1,3,1);plot(param.mu,'linewidth',2); legend('action L','action R'); xlim([1 10000]); xlabel('Trials')
ylabel('mu (environment probability estimation)');

subplot(1,3,2);plot(param.mu(:,1),param.epi(:,1),'linewidth',2);
hold on; plot(param.mu(:,2),param.epi(:,2),'linewidth',2);
xlabel('mu (environment probability estimation)'); ylabel('var/(mu(1-mu))')

subplot(1,3,3);plot(param.var(:,1),param.epi(:,1),'linewidth',2.5);
hold on; plot(param.var(:,2),param.epi(:,2),'linewidth',1.5);
xlabel('var (environment uncertainty)'); ylabel('var/(mu(1-mu))')

%%

mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT_bin40(mice);
%%
for i = 1:length(allMouse)
    [param1,param2] = upDateMouseChoice(allMouse{i}.behav);
    figure; 
    subplot(2,3,1);plot(param1.mu,'linewidth',2); legend('action L','action R'); xlabel('Trials')
    ylabel('mu (environment probability estimation)');
    
    subplot(2,3,2);plot(param1.mu(:,1),param1.epi(:,1),'linewidth',2);
    hold on; plot(param1.mu(:,2),param1.epi(:,2),'linewidth',2);
    xlabel('mu (environment probability estimation)'); ylabel('var/(mu(1-mu))')
    
    subplot(2,3,3);plot(param1.var(:,1),param1.epi(:,1),'linewidth',2.5);
    hold on; plot(param1.var(:,2),param1.epi(:,2),'linewidth',1.5);
    xlabel('var (environment uncertainty)'); ylabel('var/(mu(1-mu))')

    subplot(2,3,4);plot(param2.mu,'linewidth',2); legend('action L','action R'); xlabel('Trials')
    ylabel('mu (environment probability estimation)');
    
    subplot(2,3,5);plot(param2.mu(:,1),param2.epi(:,1),'linewidth',2);
    hold on; plot(param2.mu(:,2),param2.epi(:,2),'linewidth',2);
    xlabel('mu (environment probability estimation)'); ylabel('var/(mu(1-mu))')
    
    subplot(2,3,6);plot(param2.var(:,1),param2.epi(:,1),'linewidth',2.5);
    hold on; plot(param2.var(:,2),param2.epi(:,2),'linewidth',1.5);
    xlabel('var (environment uncertainty)'); ylabel('var/(mu(1-mu))')


end


%% all functions
function [param1,param2] = upDateMouseChoice(behav)
    startingFactor = 1; updateFactor = 0.005; 
    param.alpha = [1 1] * startingFactor; param.beta = [1 1] * startingFactor; 
    param.mu = param.alpha ./ (param.alpha + param.beta);
    param.var = param.alpha .* param.beta ./ (( param.alpha+param.beta +1) .* ( param.alpha+param.beta ).^2);
    param.epi = param.var./ (param.mu.*(1-param.mu));
    nTrials = size(behav,1); param1 = param; param2 = param;
    for t = 1:nTrials
        if behav.stimulus(t) == 1
            param1 = update(param1, behav.action(t), double(behav.action(t)==behav.stimulus(t)),t,updateFactor);
            param2 = noUpdate(param2,t);
        else
            param2 = update(param2, behav.action(t), double(behav.action(t)==behav.stimulus(t)),t,updateFactor);
            param1 = noUpdate(param1,t);
        end
    
    end



end


function [choiceP, choice, outcome] = generateChoice(param,t)
choiceP = fn_logistic(param.mu(t,2)-param.mu(t,1),4);
if choiceP > rand(); choice = 2; else; choice = 1;end
%outcome = choice-1;

if mod(t,1000)<500; outcome = choice-1;else; outcome = -choice+2; end
end

function param = update(param, action, outcome,t,updateFactor)
    if action == 1
        param.alpha(t+1,:) = [param.alpha(t,1) + outcome*updateFactor,  param.alpha(t,2)];
        param.beta(t+1,:) = [param.beta(t,1) + (1 - outcome)*updateFactor,  param.beta(t,2)];
    elseif action == 2
        param.alpha(t+1,:) = [param.alpha(t,1), param.alpha(t,2)  + outcome*updateFactor];
        param.beta(t+1,:) = [param.beta(t,1), param.beta(t,2) + (1 - outcome)*updateFactor];
    end
    param.mu(t+1,:) = param.alpha(t+1,:) ./ (param.alpha(t+1,:) + param.beta(t+1,:));
    param.var(t+1,:) = param.alpha(t+1,:) .* param.beta(t+1,:) ./ ((param.alpha(t+1,:) + param.beta(t+1,:)+1) .*...
        ( param.alpha(t+1,:)+ param.beta(t+1,:) ).^2);
    param.epi(t+1,:) = param.var(t+1,:)./ (param.mu(t+1,:).*(1-param.mu(t+1,:)));
end

function param = noUpdate(param,t)
    
    param.alpha(t+1,:) = param.alpha(t,:);
    param.beta(t+1,:) = param.beta(t,:);

    param.mu(t+1,:) = param.alpha(t+1,:) ./ (param.alpha(t+1,:) + param.beta(t+1,:));
    param.var(t+1,:) = param.alpha(t+1,:) .* param.beta(t+1,:) ./ ((param.alpha(t+1,:) + param.beta(t+1,:)+1) .*...
        ( param.alpha(t+1,:)+ param.beta(t+1,:) ).^2);
    param.epi(t+1,:) = param.var(t+1,:)./ (param.mu(t+1,:).*(1-param.mu(t+1,:)));
end