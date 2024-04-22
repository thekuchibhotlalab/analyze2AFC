function fn_runModel(state,param)
    Q = [0 0;0 0];
    histA = [0 0];
    for i = 1:size(state,1)
        [p] = fn_action_state(state,Q, histA,param);

    end



end