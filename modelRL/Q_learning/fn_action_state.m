function [p] = fn_action_state(state,Q, histA,param)
% state is 1*2 vector of [0 1]
% Q is a 2*2 vector of [L1, R1; L2, R2]
if isempty(param)
    param.alphaPos = [0.2 0.2];
    param.alphaNeg = [0.2 0.2];
    param.gammaLR = [0.2 0.2]; %size 1*2
    param.gammaRand = 0.2;
    param.gammaSwit = [-0.1 -0.1];
    param.lambdaHist = [0.9 0.9];
end

tempW = (param.gammaLR + param.gammaRand + param.gammaSwit .* histA + state * Q);

p = softmax(tempW');



end