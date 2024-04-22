function [Q] = fn_Q_update(state,Q,action,outcome,param)
 
    if outcome>0
        state_action_sel = state' * (action.* param.alphaPos);
    else
        state_action_sel = state' * (action.* param.alphaNeg);
    end

    Q = (1-state_action_sel) .* Q + state_action_sel .* outcome; 

end