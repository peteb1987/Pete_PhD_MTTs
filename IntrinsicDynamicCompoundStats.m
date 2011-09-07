function [ A_compound, Q_compound ] = IntrinsicDynamicCompoundStats( L, init_state, init_var)
%INTRINSICDYNAMICCOMPOUNDCOVARIANCE Calculate an approximation of the state
% covariance over multiple frames

Q_compound = init_var;
A_compound = eye(4);
state = init_state;

for tt = 1:L
    % Calculate Q_{tt+1} and A_{tt}
    [A, Q] = IntrinsicDynamicLinearise(state);
    
    % Update A and Q
    Q_compound = A*Q_compound*A'+Q;
    A_compound = A*A_compound;
    
    % Predict next state
    state = IntrinsicDynamicPredict(state);
end

end

